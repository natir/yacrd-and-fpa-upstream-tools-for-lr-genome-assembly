#!/usr/bin/env python3

import os
import re
import csv
import sys
import math
import argparse

import matplotlib.pyplot as plt

from collections import defaultdict

re_scrub_asm = re.compile(r"([^.]+)\.(.+)\.txt")

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", required=True)

    args = parser.parse_args(args)

    dataset2scrub2asm2metrics = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(int)
            )
        )
    )
    exp2metrics = defaultdict()
    
    for entry in os.scandir(args.directory):
        search = re_scrub_asm.search(entry.name)
        if search is None:
            continue
        
        dataset, suffix = re_scrub_asm.search(entry.name).groups()

        if "miniasm" not in suffix and "wdbtg2" not in suffix:
            scrubber = suffix
            asm = "no"
        else:
            last_point = suffix.rindex(".")
            scrubber = suffix[:last_point]
            asm = suffix[last_point+1:]
        
        with open(os.path.join(args.directory, entry.name)) as file_in:
            for row in csv.DictReader(file_in, delimiter='\t'):
                dataset2scrub2asm2metrics[dataset][scrubber][asm]["time"] = float(row["s"])
                dataset2scrub2asm2metrics[dataset][scrubber][asm]["memory"] = float(row["max_rss"])

    result = defaultdict(lambda: defaultdict(dict))

    for dataset in dataset2scrub2asm2metrics.keys():
        for scrub in dataset2scrub2asm2metrics[dataset].keys():
            if len(dataset2scrub2asm2metrics[dataset][scrub]["no"]) != 0:            
                result[dataset][scrub]["x"] = dataset2scrub2asm2metrics[dataset][scrub]["no"]["time"]
                result[dataset][scrub]["y"] = dataset2scrub2asm2metrics[dataset][scrub]["no"]["memory"]

            for asm in dataset2scrub2asm2metrics[dataset][scrub].keys():
                if asm == "no":
                    continue

                result[dataset][scrub+"_"+asm]["x"] = dataset2scrub2asm2metrics[dataset][scrub][asm]["time"] + dataset2scrub2asm2metrics[dataset][scrub]['no']["time"]
                result[dataset][scrub+"_"+asm]["y"] = max(dataset2scrub2asm2metrics[dataset][scrub]["no"]["memory"], dataset2scrub2asm2metrics[dataset][scrub][asm]["memory"])

    for dataset in result.keys():
        fig = plt.figure()
        ax = plt.axes()

        x_vals = list()
        y_vals = list()
        for analysis in result[dataset].keys():
            ax.plot(result[dataset][analysis]["x"], result[dataset][analysis]["y"], 'ro', label=analysis2label(analysis), marker=analysis2marker(analysis), color=analysis2color(analysis))
            x_vals.append(result[dataset][analysis]["x"])
            y_vals.append(result[dataset][analysis]["y"])

        min_x = math.floor(min(x_vals)-0.1*min(x_vals))
        max_x = math.ceil(max(x_vals)+0.1*max(x_vals))
        min_x = 1.1 if min_x < 1.1 else min_x
        
        min_y = math.floor(min(y_vals)-0.1*min(y_vals))
        max_y = math.ceil(max(y_vals)+0.1*max(y_vals))
        min_y = 1.1 if min_y < 1.1 else min_y

        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
        plt.xscale("log")
        plt.yscale("log")

        plt.xlabel("log of time (in second)")
        plt.ylabel("log of memory peak (in Mo)")

        #ax.legend()
        plt.savefig(dataset+"_cpu_memory.png")

        
def analysis2label(analysis):
    scrubber, *asm = analysis.split("_")
    result = ""

    if scrubber == "raw":
        result += "no_scrubb"
    elif scrubber == "dascrubber":
        result += "dascrubber"
    elif scrubber == "miniscrub.cpu":
        result += "miniscrub"
    elif scrubber == "4.4.yacrd":
        result += "yacrd"
    else:
        result += "yacrd.prec"

    if len(asm) > 0:
        asm = asm[0]
        if asm == "miniasm":
            result += "_miniasm"
        else:
            result += "_wgtdb2"
        
    return result
    
def analysis2marker(analysis):
    if "4.4.yacrd" in analysis:
        return "^"
    elif "4.4.precision.yacrd" in analysis:
        return "v"
    elif "dascrubber" in analysis:
        return "s"
    elif "miniscrub" in analysis:
        return "p"
    else:
        return "x"


def analysis2color(analysis):
    if "wdbtg2" in analysis:
        return "red"
    elif "miniasm" in analysis:
        return "green"
    else:
        return "blue"
        
    
if __name__ == "__main__":
    main(sys.argv[1:])
