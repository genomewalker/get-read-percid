import argparse
import sys
import gzip
import os
import logging
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull
from get_read_percid import __version__
import time
import json
from get_read_percid.defaults import mdmg_header, valid_ranks, filterBAM_header
from collections import defaultdict
import re
from itertools import chain

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it
def check_suffix(val, parser, var):
    if var == "--scale":
        units = ["K", "M"]
    else:
        units = ["K", "M", "G"]
    unit = val[-1]
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        if var == "--scale":
            if unit == "K":
                val = value * 1000
            elif unit == "M":
                val = value * 1000000
            elif unit == "G":
                val = value * 1000000000
            return val
        else:
            return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return value


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


def is_valid_filter(parser, arg, var, type="metaDMG"):
    if type == "metaDMG":
        header = mdmg_header
    elif type == "filterBAM":
        header = filterBAM_header
    arg = json.loads(arg)
    # check if the dictionary keys are in the mdmg header list
    for key in arg.keys():
        if key not in header:
            parser.error(
                f"argument {var}: Invalid value {key}.\n"
                f"Valid values are: {convert_list_to_str(header)}"
            )

    return arg


def is_valid_rank(parser, arg, var):
    arg = json.loads(arg)
    for key in arg.keys():
        if key not in valid_ranks.keys():
            parser.error(
                f"argument {var}: Invalid value {key}.\n"
                f"Valid values are: {convert_list_to_str(valid_ranks)}"
            )
    return arg


def get_ranks(parser, ranks, var):
    valid_ranks = [
        "domain",
        "kingdom",
        "lineage",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    ranks = ranks.split(",")
    # check if ranks are valid
    for rank in ranks:
        if rank not in valid_ranks:
            parser.error(
                f"argument {var}: Invalid value {rank}.\Rank has to be one of {convert_list_to_str(valid_ranks)}"
            )
        if rank == "all":
            ranks = valid_ranks[1:]
    return ranks


defaults = {
    "prefix": None,
    "sort_memory": "1G",
    "threads": 1,
    "chunk_size": None,
}

help_msg = {
    "bam": "The BAM file used to generate the metaDMG results",
    "prefix": "Prefix used for the output files",
    "reference_list": "A list of references that we want to extract from the BAM file",
    "rank": "Which taxonomic group and rank we want to get the reads extracted.",
    "sort_by_queryname": "Sort BAM file by queryname",
    "sort_by_coordinate": "Sort BAM file by coordinate",
    "only_damaged": "If set, only the reads damaged will be extracted",
    "sort_memory": "Set maximum memory per thread for sorting; suffix K/M/G recognized",
    "chunk_size": "Chunk size for parallel processing",
    "threads": "Number of threads",
    "debug": "Print debug messages",
    "version": "Print program version",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to extract references from BAM files and get read statistics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # parent_parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument(
        "-b",
        "--bam",
        type=lambda x: is_valid_file(parser, x, "--bam"),
        dest="bam",
        metavar="FILE",
        help=help_msg["bam"],
        required=True,
    )
    required.add_argument(
        "-r",
        "--reference-list",
        type=lambda x: is_valid_file(parser, x, "--reference-list"),
        dest="reference_list",
        metavar="FILE",
        help=help_msg["reference_list"],
    )
    optional.add_argument(
        "--prefix",
        type=str,
        default=defaults["prefix"],
        dest="prefix",
        metavar="STR",
        help=help_msg["prefix"],
    )
    optional.add_argument(
        "--sort-memory",
        type=lambda x: check_suffix(x, parser=parser, var="--sort-memory"),
        default=defaults["sort_memory"],
        dest="sort_memory",
        metavar="STR",
        help=help_msg["sort_memory"],
    )
    optional.add_argument(
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        metavar="INT",
        help=help_msg["threads"],
    )
    optional.add_argument(
        "--chunk-size",
        type=lambda x: int(
            check_values(x, minval=1, maxval=100000, parser=parser, var="--chunk-size")
        ),
        default=defaults["chunk_size"],
        dest="chunk_size",
        metavar="INT",
        help=help_msg["chunk_size"],
    )
    optional.add_argument(
        "--debug", dest="debug", action="store_true", help=help_msg["debug"]
    )
    optional.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])
    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def create_output_files(prefix, bam, taxon=None, combined=False):
    if prefix is None:
        prefix = bam.split(".")[0]
    # create output files
    if taxon is None:
        if combined:
            out_files = {
                "fastq_combined": f"{prefix}.fastq.gz",
            }
        else:
            out_files = {
                "fastq_damaged": f"{prefix}.damaged.fastq.gz",
                "fastq_nondamaged": f"{prefix}.non-damaged.fastq.gz",
                "fastq_multi": f"{prefix}.multi.fastq.gz",
                "fastq_combined": f"{prefix}.fastq.gz",
                "fastq_discarded": f"{prefix}.discarded.fastq.gz",
            }
    else:
        if combined:
            out_files = defaultdict()
            for k, v in taxon.items():
                for i in v:
                    r = splitkeep(i, "__")
                    i = re.sub("[^0-9a-zA-Z]+", "_", r[1])
                    i = re.sub("[^0-9a-zA-Z]+", "_", i)
                    out_files[f"fastq_{k}{i}"] = f"{prefix}.{k}{i}.fastq.gz"
        else:
            out_files = defaultdict()
            for k, v in taxon.items():
                for i in v:
                    r = splitkeep(i, "__")
                    i = re.sub("[^0-9a-zA-Z]+", "_", r[1])
                    out_files[
                        f"fastq_damaged_{k}{i}"
                    ] = f"{prefix}.{k}{i}.damaged.fastq.gz"
                    out_files[
                        f"fastq_nondamaged_{k}{i}"
                    ] = f"{prefix}.{k}{i}.non-damaged.fastq.gz"
                    out_files[f"fastq_multi_{k}{i}"] = f"{prefix}.{k}{i}.multi.fastq.gz"
                    out_files[f"fastq_combined_{k}{i}"] = f"{prefix}.{k}{i}.fastq.gz"
                    out_files["fastq_discarded"] = f"{prefix}.discarded.fastq.gz"
    return out_files


# From https://stackoverflow.com/a/61436083
def splitkeep(s, delimiter):
    split = s.split(delimiter)
    return [substr + delimiter for substr in split[:-1]] + [split[-1]]


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def initializer(init_data):
    global parms
    parms = init_data


# from https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.

    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize
