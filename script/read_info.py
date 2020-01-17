#!/usr/bin/env python

import os
import csv
import sys
import subprocess

import numpy
import pandas

dataset2yacrdparam = {
    "c_elegans_pb": "g800.c4.yacrd",
    "h_sapiens_chr1_ont": "g500.c4.yacrd",
    "d_melanogaster_reads_ont": "g500.c4.yacrd",
    "ERR2531093_pb": "g5000.c3.yacrd",
    "ERR2531094_pb": "g5000.c3.yacrd",
    "ERR2531095_pb": "g5000.c3.yacrd",
    "ERR2531096_pb": "g5000.c3.yacrd",
    "ERR2564864_pb": "g5000.c3.yacrd",
    "ERR2564865_pb": "g5000.c3.yacrd",
    "ERR2603033_pb": "g5000.c3.yacrd",
    "ERR2615933_pb": "g5000.c3.yacrd",
    "ERR2615934_pb": "g5000.c3.yacrd",
    "ERR2632359_pb": "g5000.c3.yacrd",
    "ERR2651535_pb": "g5000.c3.yacrd",
    "ERR2651536_pb": "g5000.c3.yacrd",
    "ERR2672424_pb": "g5000.c3.yacrd",
    "ERR2672425_pb": "g5000.c3.yacrd",
    "ERR2681660_pb": "g5000.c3.yacrd",
    "ERR2695057_pb": "g5000.c3.yacrd",
    "ERR2695058_pb": "g5000.c3.yacrd",
    "ERR3253076_pb": "g5000.c3.yacrd",
    "ERR3253077_pb": "g5000.c3.yacrd",
    "ERR3253078_pb": "g5000.c3.yacrd",
    "ERR3253079_pb": "g5000.c3.yacrd",
    "ERR3253080_pb": "g5000.c3.yacrd",
    "ERR3253081_pb": "g5000.c3.yacrd",
    "ERR3253082_pb": "g5000.c3.yacrd",
    "ERR3253083_pb": "g5000.c3.yacrd",
    "ERR3253084_pb": "g5000.c3.yacrd",
    "ERR3253085_pb": "g5000.c3.yacrd",
    "ERR3253086_pb": "g5000.c3.yacrd",
    "ERR3253087_pb": "g5000.c3.yacrd",
    "ERR3253088_pb": "g5000.c3.yacrd",
    "ERR3253089_pb": "g5000.c3.yacrd",
    "ERR3253090_pb": "g5000.c3.yacrd",
    "ERR3253091_pb": "g5000.c3.yacrd",
    "ERR3253092_pb": "g5000.c3.yacrd",
    "ERR3253093_pb": "g5000.c3.yacrd",
    "ERR3253094_pb": "g5000.c3.yacrd",
    "ERR3253095_pb": "g5000.c3.yacrd",
    "ERR3253096_pb": "g5000.c3.yacrd",
    "ERR3253097_pb": "g5000.c3.yacrd",
    "ERR3253098_pb": "g5000.c3.yacrd",
    "ERR3253099_pb": "g5000.c3.yacrd",
    "ERR3253100_pb": "g5000.c3.yacrd",
    "ERR3253102_pb": "g5000.c3.yacrd",
    "ERR3253103_pb": "g5000.c3.yacrd",
    "ERR3253104_pb": "g5000.c3.yacrd",            
    "ERR3500074_pb": "g5000.c3.yacrd",

    "SRR8494915_ont": "g500.c4.yacrd",
    "SRR8494916_ont": "g500.c4.yacrd",
    "SRR8494917_ont": "g500.c4.yacrd",
    "SRR8494918_ont": "g500.c4.yacrd",
    "SRR8494919_ont": "g500.c4.yacrd",
    "SRR8494920_ont": "g500.c4.yacrd",
    "SRR8494921_ont": "g500.c4.yacrd",
    "SRR8494922_ont": "g500.c4.yacrd",
    "SRR8494923_ont": "g500.c4.yacrd",
    "SRR8494924_ont": "g500.c4.yacrd",
    "SRR8494935_ont": "g500.c4.yacrd",
    "SRR8494936_ont": "g500.c4.yacrd",
    "SRR8494937_ont": "g500.c4.yacrd",
    "SRR8494938_ont": "g500.c4.yacrd",
    "SRR8494939_ont": "g500.c4.yacrd",
    "SRR8494940_ont": "g500.c4.yacrd",
    "SRR8494941_ont": "g500.c4.yacrd",
    "SRR8494942_ont": "g500.c4.yacrd",
    "SRR8494943_ont": "g500.c4.yacrd",
    "SRR8494944_ont": "g500.c4.yacrd",
    "SRR8494905_pb": "g800.c4.yacrd",
    "SRR8494906_pb": "g800.c4.yacrd",
    "SRR8494907_pb": "g800.c4.yacrd",
    "SRR8494908_pb": "g800.c4.yacrd",
    "SRR8494909_pb": "g800.c4.yacrd",
    "SRR8494910_pb": "g800.c4.yacrd",
    "SRR8494911_pb": "g800.c4.yacrd", 
    "SRR8494912_pb": "g800.c4.yacrd",
    "SRR8494913_pb": "g800.c4.yacrd",
    "SRR8494914_pb": "g800.c4.yacrd",
    "SRR8494925_pb": "g800.c4.yacrd",
    "SRR8494926_pb": "g800.c4.yacrd",
    "SRR8494927_pb": "g800.c4.yacrd",
    "SRR8494928_pb": "g800.c4.yacrd",
    "SRR8494929_pb": "g800.c4.yacrd",
    "SRR8494930_pb": "g800.c4.yacrd",
    "SRR8494931_pb": "g800.c4.yacrd",
    "SRR8494932_pb": "g800.c4.yacrd",
    "SRR8494933_pb": "g800.c4.yacrd",
    "SRR8494934_pb": "g800.c4.yacrd",
}

dataset2group = {
    "c_elegans_pb": "eukaryota",
    "h_sapiens_chr1_ont": "eukaryota",
    "d_melanogaster_reads_ont": "eukaryota",

    "ERR2531093_pb": "bacteria_sequel",
    "ERR2531094_pb": "bacteria_sequel",
    "ERR2531095_pb": "bacteria_sequel",
    "ERR2531096_pb": "bacteria_sequel",
    "ERR2564864_pb": "bacteria_sequel",
    "ERR2564865_pb": "bacteria_sequel",
    "ERR2603033_pb": "bacteria_sequel",
    "ERR2615933_pb": "bacteria_sequel",
    "ERR2615934_pb": "bacteria_sequel",
    "ERR2632359_pb": "bacteria_sequel",
    "ERR2651535_pb": "bacteria_sequel",
    "ERR2651536_pb": "bacteria_sequel",
    "ERR2672424_pb": "bacteria_sequel",
    "ERR2672425_pb": "bacteria_sequel",
    "ERR2681660_pb": "bacteria_sequel",
    "ERR2695057_pb": "bacteria_sequel",
    "ERR2695058_pb": "bacteria_sequel",
    "ERR3253076_pb": "bacteria_sequel",
    "ERR3253077_pb": "bacteria_sequel",
    "ERR3253078_pb": "bacteria_sequel",
    "ERR3253079_pb": "bacteria_sequel",
    "ERR3253080_pb": "bacteria_sequel",
    "ERR3253081_pb": "bacteria_sequel",
    "ERR3253082_pb": "bacteria_sequel",
    "ERR3253083_pb": "bacteria_sequel",
    "ERR3253084_pb": "bacteria_sequel",
    "ERR3253085_pb": "bacteria_sequel",
    "ERR3253086_pb": "bacteria_sequel",
    "ERR3253087_pb": "bacteria_sequel",
    "ERR3253088_pb": "bacteria_sequel",
    "ERR3253089_pb": "bacteria_sequel",
    "ERR3253090_pb": "bacteria_sequel",
    "ERR3253091_pb": "bacteria_sequel",
    "ERR3253092_pb": "bacteria_sequel",
    "ERR3253093_pb": "bacteria_sequel",
    "ERR3253094_pb": "bacteria_sequel",
    "ERR3253095_pb": "bacteria_sequel",
    "ERR3253096_pb": "bacteria_sequel",
    "ERR3253097_pb": "bacteria_sequel",
    "ERR3253098_pb": "bacteria_sequel",
    "ERR3253099_pb": "bacteria_sequel",
    "ERR3253100_pb": "bacteria_sequel",
    "ERR3253102_pb": "bacteria_sequel",
    "ERR3253103_pb": "bacteria_sequel",
    "ERR3253104_pb": "bacteria_sequel",            
    "ERR3500074_pb": "bacteria_sequel",
    
    "SRR8494915_ont": "bacteria_ont",
    "SRR8494916_ont": "bacteria_ont",
    "SRR8494917_ont": "bacteria_ont",
    "SRR8494918_ont": "bacteria_ont",
    "SRR8494919_ont": "bacteria_ont",
    "SRR8494920_ont": "bacteria_ont",
    "SRR8494921_ont": "bacteria_ont",
    "SRR8494922_ont": "bacteria_ont",
    "SRR8494923_ont": "bacteria_ont",
    "SRR8494924_ont": "bacteria_ont",
    "SRR8494935_ont": "bacteria_ont",
    "SRR8494936_ont": "bacteria_ont",
    "SRR8494937_ont": "bacteria_ont",
    "SRR8494938_ont": "bacteria_ont",
    "SRR8494939_ont": "bacteria_ont",
    "SRR8494940_ont": "bacteria_ont",
    "SRR8494941_ont": "bacteria_ont",
    "SRR8494942_ont": "bacteria_ont",
    "SRR8494943_ont": "bacteria_ont",
    "SRR8494944_ont": "bacteria_ont",
    "SRR8494905_pb": "bacteria_rsII",
    "SRR8494906_pb": "bacteria_rsII",
    "SRR8494907_pb": "bacteria_rsII",
    "SRR8494908_pb": "bacteria_rsII",
    "SRR8494909_pb": "bacteria_rsII",
    "SRR8494910_pb": "bacteria_rsII",
    "SRR8494911_pb": "bacteria_rsII", 
    "SRR8494912_pb": "bacteria_rsII",
    "SRR8494913_pb": "bacteria_rsII",
    "SRR8494914_pb": "bacteria_rsII",
    "SRR8494925_pb": "bacteria_rsII",
    "SRR8494926_pb": "bacteria_rsII",
    "SRR8494927_pb": "bacteria_rsII",
    "SRR8494928_pb": "bacteria_rsII",
    "SRR8494929_pb": "bacteria_rsII",
    "SRR8494930_pb": "bacteria_rsII",
    "SRR8494931_pb": "bacteria_rsII",
    "SRR8494932_pb": "bacteria_rsII",
    "SRR8494933_pb": "bacteria_rsII",
    "SRR8494934_pb": "bacteria_rsII",
}

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
