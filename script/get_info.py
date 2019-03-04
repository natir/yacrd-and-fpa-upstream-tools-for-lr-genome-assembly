#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse

from collections import defaultdict

import pysam

TEMPLATE_SCRUBBING = \
"""
|                | {:>30} |
| -------------- | ------------------------------:|
| \# of read     | {:>30} |
| \# of base     | {:>30} |
| N10            | {:>30} |
| L10            | {:>30} |
| N50            | {:>30} |
| L50            | {:>30} |
| N90            | {:>30} |
| L90            | {:>30} |
| % base removed | {:>30.2f} |
| \# mapped read | {:>30} |
| \# mismatch    | {:>30} |
| time           | {:>30} |
| memory         | {:>30} |
"""

TEMPLATE_CORRECTING = \
"""
|                | {:>30} |
| -------------- | ------------------------------:|
| \# of read     | {:>30} |
| \# of base     | {:>30} |
| N10            | {:>30} |
| L10            | {:>30} |
| N50            | {:>30} |
| L50            | {:>30} |
| N90            | {:>30} |
| L90            | {:>30} |
| % base removed | {:>30.2f} |
| \# mapped read | {:>30} |
| \# mismatch    | {:>30} |
| time           | {:>30} |
| memory         | {:>30} |
"""



def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--technologys", required=True, choices=['ont', 'pb'], nargs="*")
    parser.add_argument("-s", "--scrubbers", nargs="*")
    parser.add_argument("-c", "--correctors", nargs="*")
    parser.add_argument("-a", "--assemblys", nargs="*")

    args = parser.parse_args(args)

    if args.scrubbers is None and args.correctors is None and args.assemblys is None:
        print("Noting to do", file=sys.stderr)
    elif args.correctors is None and args.assemblys is None:
        data = stat_scrubber(args.technologys, args.scrubbers)
        show_scrubber(data)
    elif args.assemblys is None:
        data = stat_correction(args.technologys, args.scrubbers, args.correctors)
        show_scrubber(data)
    else:
        stat_correction(args.technologys, args.scrubbers, args.correctors, args.assemblys)


def stat_scrubber(technologys, scrubbers):
    data = {"name": list(), "# read": list(), "# base": list(), "coverage": list(),"N10": list(), "N50": list(), "N90": list(), "L10": list(), "L50": list(), "L90": list(), "% base removed": list(), "# mapped read": list(), "# mismatch": list(), "time": list(), "memory": list()}
    
    for tech in technologys:
        for scrub in scrubbers:
            __stat_scrubber(tech, scrub, data)

    return data

def show_scrubber(data):
    print("|                | " + " | ".join(data["name"]) + " |")
    print("| - |"+" -:| "*len(data["name"]))
    print("| \# of read     | " + " | ".join(data["# read"]) + " |")
    print("| \# of base     | " + " | ".join(data["# base"]) + " |")
    print("| N10            | " + " | ".join(data["N10"]) + " |")
    print("| L10            | " + " | ".join(data["L10"]) + " |")
    print("| N50            | " + " | ".join(data["N50"]) + " |")
    print("| L50            | " + " | ".join(data["N50"]) + " |")
    print("| N90            | " + " | ".join(data["N90"]) + " |")
    print("| L90            | " + " | ".join(data["L90"]) + " |")
    print("| % base removed | " + " | ".join(data["% base removed"]) + " |")
    print("| \# mapped read | " + " | ".join(data["# mapped read"]) + " |")
    print("| \# mismatch    | " + " | ".join(data["# mismatch"]) + " |")
    print("| time           | " + " | ".join(data["time"]) + " |")
    print("| memory         | " + " | ".join(data["memory"]) + " |")
    
def __stat_scrubber(tech, scrub, data):
    raw_reads      = "data/real_reads_{}.fasta"
    benchmark_file = "benchmarks/real_reads_{}.{}.txt"
    scrubbed_reads = "scrubbing/real_reads_{}.{}.fasta"
    mapping_file   = "mapping/scrubbing/real_reads_{}.{}.bam"

    s_nb_base, s_lengths = reads_stat(scrubbed_reads.format(tech, scrub))

    raw_nb_base, _ = reads_stat(raw_reads.format(tech))

    percent_base = 100 - (s_nb_base * 100 / raw_nb_base)

    s_lengths.reverse()
    n10, l10 = get_N_L(s_lengths, 0.1)
    n50, l50 = get_N_L(s_lengths, 0.5)
    n90, l90 = get_N_L(s_lengths, 0.9)

    nb_map, edit_distance_sum = parse_mapping(mapping_file.format(tech, scrub))
    
    
    if os.path.isfile(benchmark_file.format(tech, scrub)):
        time, memory = get_benchmark(benchmark_file.format(tech, scrub))
    else:
        time, memory = 0, 0

    data["name"].append(generate_column_name(tech, scrub))
    data["# read"].append(str(len(s_lengths)))
    data["# base"].append(str(s_nb_base))
    data["coverage"].append("{:.2f}x".format(s_nb_base / 5231428))
    data["N10"].append(str(n10))
    data["N50"].append(str(n50))
    data["N90"].append(str(n90))
    data["L10"].append(str(l10))
    data["L50"].append(str(l50))
    data["L90"].append(str(l90))
    data["% base removed"].append("{:.2f}".format(percent_base))
    data["# mapped read"].append(str(nb_map))
    data["# mismatch"].append(str(edit_distance_sum))
    data["time"].append(str(time))
    data["memory"].append(str(memory))

    
def stat_correction(techs, scrubs, corrs):
    data = {"name": list(), "# read": list(), "# base": list(), "N10": list(), "N50": list(), "N90": list(), "L10": list(), "L50": list(), "L90": list(), "% base removed": list(), "# mapped read": list(), "# mismatch": list(), "time": list(), "memory": list()}
    

    for tech in techs:
        for scrub in scrubs:
            for corr in corrs:
                __stat_correction(tech, scrub, corr, data)

    return data


def __stat_correction(tech, scrub, corr, data):
    raw_reads      = "scrubbing/real_reads_{}.{}.fasta"
    benchmark_file = "benchmarks/real_reads_{}.{}.{}.txt"
    corrected_reads = "correction/real_reads_{}.{}.{}.fasta"
    mapping_file   = "mapping/correction/real_reads_{}.{}.{}.bam"

    if not os.path.isfile(corrected_reads.format(tech, scrub, corr)):
        return 
    
    s_nb_base, s_lengths = reads_stat(corrected_reads.format(tech, scrub, corr))

    raw_nb_base, _ = reads_stat(raw_reads.format(tech, scrub))

    percent_base = 100 - (s_nb_base * 100 / raw_nb_base)

    s_lengths.reverse()
    n10, l10 = get_N_L(s_lengths, 0.1)
    n50, l50 = get_N_L(s_lengths, 0.5)
    n90, l90 = get_N_L(s_lengths, 0.9)
    
    if os.path.isfile(benchmark_file.format(tech, scrub, corr)):
        time, memory = get_benchmark(benchmark_file.format(tech, scrub, corr))
    else:
        time, memory = 0, 0

    nb_map, edit_distance_sum = parse_mapping(mapping_file.format(tech, scrub, corr))
    
    data["name"].append(generate_column_name(tech, scrub, corr))
    data["# read"].append(str(len(s_lengths)))
    data["# base"].append(str(s_nb_base))
    data["N10"].append(str(n10))
    data["N50"].append(str(n50))
    data["N90"].append(str(n90))
    data["L10"].append(str(l10))
    data["L50"].append(str(l50))
    data["L90"].append(str(l90))
    data["% base removed"].append("{:.2f}".format(percent_base))
    data["# mapped read"].append(str(nb_map))
    data["# mismatch"].append(str(edit_distance_sum))
    data["time"].append(str(time))
    data["memory"].append(str(memory))


def stat_assembly(technology, scrubber, corrector, assembly):
    pass        


def get_N_L(lengths, part):
    all_base = sum(lengths)

    index = -1
    counts = 0
    for i, val in enumerate(lengths):
        if counts >= (all_base * part):
            index = i
            break
        
        counts += val

    if index == -1:
        return 0, 0
    
    return index, lengths[index]

def parse_mapping(mapping_file):
    edit_distance = 0
    read_id = set()

    mapping = pysam.AlignmentFile(mapping_file, "rb")
    for m in mapping.fetch():
        if m.flag == 0 or m.flag == 16:
            read_id.add(m.query_name)
            edit_distance += m.get_tag("NM")
        
    return len(read_id), edit_distance
        

def get_benchmark(benchmark_file):
    with open(benchmark_file) as benchmark_handler:
        reader = csv.DictReader(benchmark_handler, delimiter='\t')
        for row in reader:
            return row["s"], row["max_rss"]
        
            
def generate_column_name(*name):
    return ".".join(name)

def reads_stat(path):
    lengths = list()
    nb_base = 0
    
    with open(path) as file_handler:
        for line in file_handler:
            if line.startswith(">"):
                continue

            length = len(line.strip())
            
            nb_base += length
            lengths.append(length)

    lengths.sort()
            
    return nb_base, lengths


if __name__ == "__main__":
    main(sys.argv[1:])
