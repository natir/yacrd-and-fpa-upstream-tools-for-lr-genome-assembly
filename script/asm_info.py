#!/usr/bin/env python

import os
import csv
import sys
import subprocess

import numpy
import pandas

from utils_for_info import *

def main(args=None):
    
    if args is None:
        args = sys.argv[1:]

    index = pandas.MultiIndex(levels=[[],[], [], []], codes=[[], [], [], []], names=["dataset", "scrubber", "assembly", "group"])

    columns = [
               "#contigs",
               "NGA50",
               "NGA50 10kb",
               "Largest contig",
               "Largest alignment",
               "Total length",
               "Reference length",
               "Asm/Ref",
               "Indels per 100kb",
               "Mismatches per 100kb",
               "#relocations",
               "#translocations",
               "#inversions",
               "cumulative len of relocations",
        ]

    df = pandas.DataFrame(index=index, columns=columns)

    for dataset in dataset2group.keys():
        for scrubber in ["raw", "dascrubber", "yacrd"]:
            for assembly in ["miniasm", "wtdbg2"]:
                if scrubber == "yacrd":
                    scrubber_param = dataset2yacrdparam[dataset]
                else:
                    scrubber_param = scrubber

                if scrubber == "dascrubber" and dataset in dascrubber_skip:
                    continue

                if dataset in quast_skip:
                    continue
                    
                default = get_asm_info("quast", dataset, scrubber_param, assembly)
                large_mis_assembly = get_asm_info("quast_mis_size_10000", dataset, scrubber_param, assembly)
                mis_assembly = get_mis_assembly_info("quast_mis_size_10000", dataset, scrubber_param, assembly)

                df.loc[(clean_name(dataset), scrubber, assembly, dataset2group[dataset]),:] = [
                                                                                               default["# contigs"],
                                                                                               default["NGA50"],
                                                                                               large_mis_assembly["NGA50"],
                                                                                               default["Largest contig"],
                                                                                               default["Largest alignment"],
                                                                                               default["Total length"],
                                                                                               default["Reference length"],
                                                                                               default["Total length"] / default["Reference length"],
                                                                                               default["# indels per 100 kbp"],
                                                                                               default["# mismatches per 100 kbp"],
                                                                                               mis_assembly["# c. relocations"],
                                                                                               mis_assembly["# c. translocations"],
                                                                                               mis_assembly["# c. inversions"],
                                                                                               get_cumulative_len_of_relocation("quast", dataset, scrubber_param, assembly),
                ]
                
    group_dataset_asm = df.groupby(level=[0, 2])
    df_ratio = pandas.DataFrame(index=index, columns=columns)

    for key, item in group_dataset_asm:
        group = item.index.values[0][3]

        if (key[0], "dascrubber", key[1], group) in item.index:
            df_ratio.loc[(key[0], "dascrubber", key[1], group),:] = (item.loc[(key[0], "dascrubber",)] / item.loc[(key[0], "raw",)]).loc[(key[1], group),].tolist()
        
        if (key[0], "yacrd", key[1], group) in item.index:
            df_ratio.loc[(key[0], "yacrd", key[1], group),:] = (item.loc[(key[0], "yacrd",)] / item.loc[(key[0], "raw",)]).loc[(key[1], group),].tolist()


    if len(args) == 1 and args[0] == "latex":
        print(df.reset_index(level=3, drop=True).to_latex())
        print(df_ratio.reset_index(level=3, drop=True).to_latex())
        print(df.groupby(level=[1, 2, 3]).mean().to_latex())
    else:
        print(df.reset_index(level=3, drop=True).to_csv())
        print(df_ratio.reset_index(level=3, drop=True).to_csv())
        print(df.groupby(level=[1, 2, 3]).mean().to_csv())
        
    
def clean_name(dataset_name):
    return "_".join(dataset_name.split("_")[:-1])
        
def get_asm_info(prefix, dataset, scrubber, assembly):
    return get_quast_data(f"{prefix}/{dataset}.{scrubber}.{assembly}/report.tsv")
    
def get_quast_data(filename):
    quast_data = dict()

    quast_data = {
        "# contigs": numpy.nan,
        "Largest contig": numpy.nan,
        "Total length": numpy.nan,
        "Reference length": numpy.nan,
        "NGA50": numpy.nan,
        "Largest alignment": numpy.nan,
        "# mismatches per 100 kbp": numpy.nan,
        "# indels per 100 kbp": numpy.nan,
        "# misassemblies": numpy.nan, 
    }
    
    try:
        with open(filename) as file_handler:
            reader = csv.reader(file_handler, delimiter='\t')
            for row in reader:
                if row[0] in quast_data:
                    try:
                        quast_data[row[0]] = float(row[1])
                    except ValueError:
                        quast_data[row[0]] = numpy.nan
    except FileNotFoundError:
        print(f"Error can't open quast report file {filename}", file=sys.stderr)

    return quast_data

def get_mis_assembly_info(prefix, dataset, scrubber, assembly):
    filename = f"{prefix}/{dataset}.{scrubber}.{assembly}/contigs_reports/misassemblies_report.tsv"

    data = {
        "# c. relocations": numpy.nan,
        "# c. translocations": numpy.nan,
        "# c. inversions": numpy.nan,
    }

    try:
        with open(filename) as file_handler:
            reader = csv.reader(file_handler, delimiter='\t')

            for row in reader:
                if row[0].strip() in data:
                    data[row[0].strip()] = float(row[1])
    except FileNotFoundError:
        print(f"Error can't open mis assembly file {filename}", file=sys.stderr)

    return data

def get_cumulative_len_of_relocation(prefix, dataset, scrubber, assembly):
    filename = f"{prefix}/{dataset}.{scrubber}.{assembly}/contigs_reports/all_alignments_{dataset}-{scrubber.replace('.', '-')}-{assembly}.tsv"

    try:
        cum_relocation_len = 0
        with open(filename) as file_handler:
            for line in file_handler:
                if line.startswith("relocation"):
                    relocation_len = line.split("=")[1].strip().split(" ")[0].strip()
                    cum_relocation_len += abs(int(relocation_len))

            return cum_relocation_len
    except FileNotFoundError:
        print(f"Error can't open misassembly report {filename}", file=sys.stderr)
        return numpy.nan

if __name__ == "__main__":
    main(sys.argv[1:])
