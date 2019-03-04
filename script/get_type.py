#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse

from collections import defaultdict

raw_path = "data/real_reads_{}.fasta"
yacrd_path = "scrubbing/real_reads_{}.yacrd2"
miniscrub_path = "scrubbing/real_reads_{}.miniscrub.fasta"
dascrubber_path = "scrubbing/real_reads_{}.dascrubber.fasta"
def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-t", "--technology", required=True, choices=['ont', 'pb'])
    
    args = parser.parse_args(args)

    y_discard, y_splited, y_trimmed, y_nmodified = yacrd_analysis(yacrd_path.format(args.technology), raw_path.format(args.technology))
    m_discard, m_splited, m_trimmed, m_nmodified = miniscrub_analysis(miniscrub_path.format(args.technology), raw_path.format(args.technology))
    d_discard, d_splited, d_trimmed, d_nmodified = dascrubber_analysis(dascrubber_path.format(args.technology), raw_path.format(args.technology))

    print("|           | yacrd | dascrubber | miniscrub |")
    print("| --------- | -----:| ----------:| ---------:|")
    print("| discard   | {:5} | {:10} | {:9} |".format(len(y_discard), len(m_discard), len(d_discard)))
    print("| splited   | {:5} | {:10} | {:9} |".format(len(y_splited), len(m_splited), len(d_splited)))
    print("| trimmed   | {:5} | {:10} | {:9} |".format(len(y_trimmed), len(m_trimmed), len(d_trimmed)))
    print("| nmodified | {:5} | {:10} | {:9} |".format(len(y_nmodified), len(m_nmodified), len(d_nmodified)))


    print("\ndiscard ")
    print("|            | yacrd | dascrubber | miniscrub |")
    print("| ---------- | -----:| ----------:| ---------:|")
    print("| yacrd      |       |            |           |")
    print("| dascrubber | {:5.2f} |            |           |".format(len(y_discard & d_discard) / len(y_discard | d_discard)))
    print("| miniscrub  | {:5.2f} | {:10.2f} |           |".format(len(y_discard & m_discard) / len(y_discard | m_discard), len(m_discard & d_discard) / len(m_discard | d_discard)))

    print("\nsplited")
    print("|            | yacrd | dascrubber | miniscrub |")
    print("| ---------- | -----:| ----------:| ---------:|")
    print("| yacrd      |       |            |           |")
    print("| dascrubber | {:5.2f} |            |           |".format(len(y_splited & d_splited) / len(y_splited | d_splited)))
    print("| miniscrub  | {:5.2f} | {:10.2f} |           |".format(len(y_splited & m_splited) / len(y_splited | m_splited), len(m_splited & d_splited) / len(m_splited | d_splited)))

    print("\ntrimmed")
    print("|            | yacrd | dascrubber | miniscrub |")
    print("| ---------- | -----:| ----------:| ---------:|")
    print("| yacrd      |       |            |           |")
    print("| dascrubber | {:5.2f} |            |           |".format(len(y_trimmed & d_trimmed) / len(y_trimmed | d_trimmed)))
    print("| miniscrub  | {:5.2f} | {:10.2f} |           |".format(len(y_trimmed & m_trimmed) / len(y_trimmed | m_trimmed), len(m_trimmed & d_trimmed) / len(m_trimmed | d_trimmed)))

    print("\nnot modified")
    print("|            | yacrd | dascrubber | miniscrub |")
    print("| ---------- | -----:| ----------:| ---------:|")
    print("| yacrd      |       |            |           |")
    print("| dascrubber | {:5.2f} |            |           |".format(len(y_nmodified & d_nmodified) / len(y_nmodified | d_nmodified)))
    print("| miniscrub  | {:5.2f} | {:10.2f} |           |".format(len(y_nmodified & m_nmodified) / len(y_nmodified | m_nmodified), len(m_nmodified & d_nmodified) / len(m_nmodified | d_nmodified)))
    
    
def yacrd_analysis(filepath, rawpath):
    discard, splited, trimmed, nmodified = set(), set(), set(), set()

    with open(filepath) as yacrd_report:
        reader = csv.reader(yacrd_report, delimiter='\t')
        for row in reader:
            if row[0] == "Not_covered":
                discard.add(row[1])
            elif row[0] == "Chimeric":
                splited.add(row[1])
            elif row[0] == "NotBad":
                trimmed.add(row[1])
            else:
                pass
                #print(row)
                
    with open(rawpath) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                ident = line.strip().split()[0][1:]
                if ident not in discard and ident not in splited and ident not in trimmed:
                    nmodified.add(ident)
                
    return discard, splited, trimmed, nmodified


def miniscrub_analysis(filepath, rawpath):
    discard, splited, trimmed, nmodified = set(), set(), set(), set()

    reads2length = dict()
    ident = ""
    with open(rawpath) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
              ident = line.strip().split()[0][1:]
            else:
                reads2length[ident] = len(line.strip())

    reads2cuts = defaultdict(list)
    discard = set(reads2length.keys())
    with open(filepath) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                ident = line.strip().split()[0][1:]
                if "_" in ident:
                    ident, cut = ident.split("_")
                    discard.discard(ident)
                    reads2cuts[ident].append(cut)
                else:
                    discard.discard(ident)
                    nmodified.add(ident)

    for k, v in reads2cuts.items():
        if len(v) > 1:
            splited.add(k)
        else:
            trimmed.add(k)

    return discard, splited, trimmed, nmodified


def dascrubber_analysis(filepath, rawpath):
    discard, splited, trimmed, nmodified = set(), set(), set(), set()

    reads2length = dict()
    ident = ""
    with open(rawpath) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
              ident = line.strip().split()[0][1:]
            else:
                reads2length[ident] = len(line.strip())

    reads2cuts = defaultdict(list)
    discard = set(reads2length.keys())
    with open(filepath) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                ident = line.strip().split()[0][1:]
                if "_" in ident:
                    ident, cut = ident.split("/")
                    discard.discard(ident)
                    reads2cuts[ident].append(cut)
               
    for k, v in reads2cuts.items():
        if len(v) > 1:
            splited.add(k)
        else:
            begin, end = v[0].split("_")
            if begin == "0" and end == str(reads2length[k]):
                nmodified.add(k)
            trimmed.add(k)

    return discard, splited, trimmed, nmodified


if __name__ == "__main__":
    main(sys.argv[1:])
