import tqdm
import pysam
import logging
from multiprocessing import Pool
from functools import partial
from get_read_percid.utils import is_debug, calc_chunksize, initializer
import numpy as np
import sys
import gzip

# import cProfile as profile
# import pstats

log = logging.getLogger("my_logger")


def get_percid(bam, prefix=None, threads=1):
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    # write results to csv file
    if prefix is not None:
        ofname = f"{prefix}.percid.csv.gz"
    else:
        ofname = "sample.percid.csv.gz"
    # get number of reads in the file
    n_reads = samfile.count(until_eof=True)
    samfile.close()
    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)
    _open = partial(gzip.open, mode="wt") if ofname.endswith(".gz") else open
    with _open(ofname) as f:
        for aln in tqdm.tqdm(
            samfile.fetch(until_eof=True),
            total=n_reads,
            desc="Processing",
            disable=is_debug(),
            unit=" reads",
        ):
            ani_read = (1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100
            f.write(f"{aln.query_name},{aln.reference_name},{ani_read}\n")
    samfile.close()


def get_ids_and_filter(params, refs, reads, prefix, threads=1):
    bam, references = params

    samfile = pysam.AlignmentFile(bam, "rb", threads=threads)

    # convert the dictionary to an array
    refs_dict = dict(zip(samfile.references, samfile.lengths))
    for reference in references:
        my_array = np.array(list(refs_dict.items()))

        # get the indices of the keys to keep
        keep_indices = np.isin(my_array[:, 0], [reference])

        # use array indexing to get the key-value pairs to keep
        filtered_array = my_array[keep_indices]

        # convert the filtered array back to a dictionary
        refs_dict = dict(filtered_array)
        refs_dict = {refs[x]: y for x, y in refs_dict.items()}

        (ref_names, ref_lengths) = zip(*refs_dict.items())
        ref_lengths = list(ref_lengths)
        # convert reference lengths to integers
        ref_lengths = [int(x) for x in ref_lengths]

        refs_idx = {sys.intern(str(x)): i for i, x in enumerate(ref_names)}

        ref_name = refs[reference]
        # write results to csv file
        if prefix is not None:
            ofname = f"{prefix}___{ref_name}--percid.csv.gz"
            obam_file = f"{prefix}___{ref_name}.bam"
        else:
            ofname = f"{ref_name}--percid.csv.gz"
            obam_file = f"{ref_name}.bam"

        out_bam_file = pysam.AlignmentFile(
            obam_file,
            "wb",
            referencenames=list(ref_names),
            referencelengths=ref_lengths,
            threads=1,
        )

        _open = partial(gzip.open, mode="wt") if ofname.endswith(".gz") else open
        with _open(ofname) as f:
            for aln in samfile.fetch(
                contig=reference, multiple_iterators=False, until_eof=True
            ):
                if reads is None:
                    ani_read = (
                        1 - ((aln.get_tag("NM") / aln.infer_query_length()))
                    ) * 100
                    f.write(f"{aln.query_name},{ani_read}\n")
                    aln.reference_id = refs_idx[ref_name]
                    out_bam_file.write(aln)
                else:
                    if aln.query_name in reads:
                        if reads[aln.query_name] == reference:
                            ani_read = (
                                1 - ((aln.get_tag("NM") / aln.infer_query_length()))
                            ) * 100
                            f.write(f"{aln.query_name},{ani_read}\n")
                            aln.reference_id = refs_idx[ref_name]
                            out_bam_file.write(aln)
        out_bam_file.close()

    samfile.close()


def get_read_percid(
    bam,
    refs,
    reads,
    prefix=None,
    chunksize=None,
    threads=1,
):
    # prof = profile.Profile()
    # prof.enable()
    ref_list = [key for key in refs]

    if (chunksize is not None) and ((len(ref_list) // chunksize) > threads):
        c_size = chunksize
    else:
        c_size = calc_chunksize(n_workers=threads, len_iterable=len(ref_list), factor=4)

    ref_chunks = [ref_list[i : i + c_size] for i in range(0, len(ref_list), c_size)]

    params = zip([bam] * len(ref_chunks), ref_chunks)

    if is_debug():
        data = list(
            map(
                partial(
                    get_ids_and_filter,
                    refs=refs,
                    reads=reads,
                    prefix=prefix,
                    threads=threads,
                ),
                params,
            )
        )
    else:
        logging.info(
            f"Processing {len(ref_chunks):,} chunks of {c_size:,} references each..."
        )
        p = Pool(
            threads,
            initializer=initializer,
            initargs=([params, refs],),
        )

        data = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(
                        get_ids_and_filter,
                        refs=refs,
                        reads=reads,
                        prefix=prefix,
                        threads=threads,
                    ),
                    params,
                    chunksize=1,
                ),
                total=len(ref_chunks),
                leave=False,
                ncols=80,
                desc="References processed",
            )
        )

        p.close()
        p.join()
