"""
 Copyright (c) 2022 Antonio Fernandez-Guerra

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 """


import logging

from get_read_percid.utils import (
    get_arguments,
    create_output_files,
    splitkeep,
    fast_flatten,
)
from get_read_percid.lib import (
    load_references,
    filter_damaged_taxa,
    load_fb_results,
    filter_fb_references,
)
from get_read_percid.defaults import valid_ranks
from get_read_percid.extract import get_read_percid, get_percid
import pandas as pd
import os
import pysam
import numpy as np
from Bio import SeqIO, SeqRecord
import gzip
from mimetypes import guess_type
from functools import partial
import tqdm
from collections import defaultdict
import re


log = logging.getLogger("my_logger")


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = get_arguments()
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )

    # Check that rank and taxonomy file are both present
    if args.reference_list:
        log.info("Loading references list...")
        references = load_references(args.reference_list)
        if references is None:
            log.error(
                "The reference list file has an invalid number of columns. Exiting"
            )
    else:
        references = None

    logging.info("Loading BAM file...")
    save = pysam.set_verbosity(0)
    bam = args.bam
    samfile = pysam.AlignmentFile(bam, "rb", threads=args.threads)

    chr_lengths = []
    for chrom in samfile.references:
        chr_lengths.append(samfile.get_reference_length(chrom))
    max_chr_length = np.max(chr_lengths)

    sorted_bam = bam.replace(".bam", ".id-sorted.bam")
    sorted_bam_index_bai = sorted_bam.replace(".id-sorted.bam", ".id-sorted.bam.bai")
    sorted_bam_index_csi = sorted_bam.replace(".id-sorted.bam", ".id-sorted.bam.csi")
    # Check if BAM files is not sorted by coordinates, sort it by coordinates
    if not samfile.header["HD"]["SO"] == "coordinate":
        log.info("BAM file is not sorted by coordinates, sorting it...")
        pysam.sort(
            "-@", str(args.threads), "-m", str(args.sort_memory), "-o", sorted_bam, bam
        )
        bam = sorted_bam
        samfile = pysam.AlignmentFile(bam, "rb", threads=args.threads)

    if samfile.has_index():
        if os.path.exists(sorted_bam_index_bai):
            os.remove(sorted_bam_index_bai)
        elif os.path.exists(sorted_bam_index_csi):
            os.remove(sorted_bam_index_csi)

    logging.info("Indexing BAM file...")
    if max_chr_length > 536870912:
        logging.info("A reference is longer than 2^29, indexing with csi")
        pysam.index("-c", "-@", str(args.threads), bam)
    else:
        pysam.index(
            "-@", str(args.threads), bam
        )  # Need to reload the samfile after creating index
        log.info("Re-loading BAM file")
        samfile = pysam.AlignmentFile(bam, "rb", threads=args.threads)
    pysam.set_verbosity(save)

    refs_bam = [
        chrom.contig for chrom in samfile.get_index_statistics() if chrom.mapped > 0
    ]
    samfile.close()

    # If we have a taxonomy file, get the references that are in the BAM file
    if references is not None:
        logging.info("Extracting stats and reads for references in list...")
        # Check if references are in the BAM file using set intersection
        refs_tax = set(references.keys()).intersection(set(refs_bam))
        refs_tax = list(refs_tax)
        # check that references in refs_tax are the same than in references
        if not set(refs_tax) == set(references.keys()):
            log.error(
                "Some references in the taxonomy file are not in the BAM file. Exiting"
            )
            exit(1)
        if not refs_tax:
            log.error("No references found in BAM file")
            exit(1)
        reads = get_read_percid(
            bam=bam,
            refs=references,
            threads=args.threads,
            prefix=args.prefix,
            chunksize=args.chunk_size,
        )
    else:
        logging.info("Extracting stats and reads for all references...")
        reads = get_percid(
            bam=bam,
            threads=args.threads,
            prefix=args.prefix,
        )

    if os.path.exists(sorted_bam_index_bai):
        os.remove(sorted_bam_index_bai)
    elif os.path.exists(sorted_bam_index_csi):
        os.remove(sorted_bam_index_csi)
    elif os.path.exists(sorted_bam):
        os.remove(sorted_bam)
    logging.info("Done!")


if __name__ == "__main__":
    main()
