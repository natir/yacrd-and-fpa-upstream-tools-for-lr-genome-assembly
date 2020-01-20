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

    index = pandas.MultiIndex(levels=[[],[], []], codes=[[],[], []], names=["dataset", "scrubber", "group"])
        
    df = pandas.DataFrame(index=index, columns=["# reads mapped", "error rate", "total length", "# chimera"])
    
    for dataset in dataset2group.keys():
        for scrubber in ["raw", "dascrubber", "yacrd"]:
            if scrubber == "yacrd":
                scrubber_param = dataset2yacrdparam[dataset]
            else:
                scrubber_param = scrubber
                
            map_info = get_mapping_info(dataset, scrubber_param)
            df.loc[(clean_name(dataset), scrubber, dataset2group[dataset]),:] = [map_info["reads mapped:"], map_info["error rate:"], map_info["total length:"], get_chimera_info(dataset, scrubber_param)]

    if len(args) == 1 and args[0] == "latex":
        print(df.reset_index(level=2, drop=True).to_latex())
        print(df.groupby(level=[1, 2]).mean().to_latex())
    else:
        print(df.reset_index(level=2, drop=True).to_csv())
        print(df.groupby(level=[1, 2]).mean().to_csv())

def clean_name(dataset_name):
    return "_".join(dataset_name.split("_")[:-1])
        
def get_chimera_info(dataset_name, scrubber):
    
    p = subprocess.Popen(["./script/found_chimera.py", f"mapping/{dataset_name}.{scrubber}.paf"], stdout=subprocess.PIPE, universal_newlines=True)

    stdout, stderr = p.communicate()

    status = p.wait()

    if status == 0:
        return int(stdout)
    else:
        print(f"Error with chimeric detection for {dataset_name}.{scrubber} return code {status}", file=sys.stderr)
        print(f"{stderr}", file=sys.stderr)
        return numpy.nan
        
def get_mapping_info(dataset_name, scrubber):
    
    p = subprocess.Popen(["samtools", "stats", f"mapping/{dataset_name}.{scrubber}.bam"], stdout=subprocess.PIPE, universal_newlines=True)

    stdout, stderr = p.communicate()

    status = p.wait()

    keept_value = {
        "reads mapped:": numpy.nan,
        "error rate:": numpy.nan,
        "total length:": numpy.nan
    }

    if status == 0:
        for line in stdout.split("\n"):
            if line.startswith("SN"):
                row = line.split("\t")
                if row[1] in keept_value:
                    keept_value[row[1]] = float(row[2])
    else:
        print(f"Error with mapping stat extraction for {dataset_name}.{scrubber} return code {status}", file=sys.stderr)
        print(f"{stderr}", file=sys.stderr)

    return keept_value
        
if __name__ == "__main__":
    main(sys.argv[1:])
