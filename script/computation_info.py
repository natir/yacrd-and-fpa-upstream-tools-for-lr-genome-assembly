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

    index = pandas.MultiIndex(levels=[[],[], [], []], codes=[[],[], [], []], names=["dataset", "scrubber", "assembly", "group"])

    df = pandas.DataFrame(index=index, columns=["time (s)", "memory (MB)", "minimap time (s)", "minimap memory (MB)", "yacrd time (s)", "yacrd memory (MB)"])

    for dataset in dataset2group.keys():
        for scrubber in ["raw", "dascrubber", "yacrd"]:
            if scrubber == "yacrd":
                minimap_param = dataset2yacrdparam[dataset].split(".")[0]
                minimap_values = get_values_from_file(f"benchmarks/{dataset}.{minimap_param}.minimap.yacrd.txt")
                yacrd_values = get_values_from_file(f"benchmarks/{dataset}.{dataset2yacrdparam[dataset]}.txt")
                
                df.loc[(clean_name(dataset), scrubber, "no", dataset2group[dataset]),:] = [
                    minimap_values["time"] + yacrd_values["time"],
                    minimap_values["mem"] + yacrd_values["mem"],
                    minimap_values["time"],
                    minimap_values["mem"],
                    yacrd_values["time"],
                    yacrd_values["mem"],
                ]
            elif scrubber == "dascrubber":
                values = get_values(dataset, scrubber)
                df.loc[(clean_name(dataset), scrubber, "no", dataset2group[dataset]),:] = [
                    values["time"],
                    values["mem"],
                    numpy.nan,
                    numpy.nan,
                    numpy.nan,
                    numpy.nan,
                ]

            for assembly in ["miniasm", "wdbtg2"]:
                if scrubber == "yacrd":
                    scrubber_param = dataset2yacrdparam[dataset]
                else:
                    scrubber_param = scrubber
                      
                values = get_values(dataset, scrubber_param, assembly)
                df.loc[(clean_name(dataset), scrubber, assembly, dataset2group[dataset]),:] = [
                    values["time"],
                    values["mem"],
                    numpy.nan,
                    numpy.nan,
                    numpy.nan,
                    numpy.nan,
                ]

    if len(args) == 1 and args[0] == "latex":
        print(df.reset_index(level=3, drop=True).to_latex())
        print(df.groupby(level=[1, 2, 3]).mean().to_latex())
    else:
        print(df.reset_index(level=3, drop=True).to_csv())
        print(df.groupby(level=[1, 2, 3]).mean().to_csv())

                
def clean_name(dataset_name):
    return "_".join(dataset_name.split("_")[:-1])

def get_values(dataset, scrubber, assembly=None):
    if assembly is None:
        filename = f"benchmarks/{dataset}.{scrubber}.txt"
    else:
        filename = f"benchmarks/{dataset}.{scrubber}.{assembly}.txt"

    return get_values_from_file(filename)

def get_values_from_file(filename):
    data = {
        "time": numpy.nan,
        "mem": numpy.nan,
    }
        
    try:
        with open(filename) as file_handler:
            reader = csv.reader(file_handler, delimiter='\t')
            next(reader) # skip first line
            row = next(reader)
            data["time"] = float(row[0])
            data["mem"] = float(row[2])
    except FileNotFoundError:
        print(f"Error can't open quast report file {filename}", file=sys.stderr)

    return data
            
if __name__ == "__main__":
    main(sys.argv[1:])
