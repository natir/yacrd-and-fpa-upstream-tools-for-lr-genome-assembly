#!/usr/bin/env python

import os
import csv
import sys
import subprocess

import numpy
import pandas

from utils_for_info import *
from asm_info import get_cumulative_len_of_relocation_from_filename


def main(args=None):
    
    if args is None:
        args = sys.argv[1:]

    index = pandas.MultiIndex(levels=[[], []], codes=[[], []], names=["dataset", "fpa"])

    columns = [
               "cumulative len of relocations",
    ]

    df = pandas.DataFrame(index=index, columns=columns)

    datasets = [
                "c_elegans_pb",
                "h_sapiens_chr1_ont",
                "d_melanogaster_reads_ont",
                "real_reads_pb",
                "real_reads_ont",
        ]

    for dataset in datasets:
        for prefix in ["", "fpa_"]:
            filename = f"fpa/quast/{prefix}{dataset}/contigs_reports/all_alignments_{prefix}{dataset}.tsv"
            df.loc[(dataset, prefix == "fpa_"),:] = [
                                                     get_cumulative_len_of_relocation_from_filename(filename)
                                                     ]

    print(df.to_csv())
            
if __name__ == "__main__":
    main(sys.argv[1:])
