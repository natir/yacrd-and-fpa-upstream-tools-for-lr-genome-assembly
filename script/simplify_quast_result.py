#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse

from collections import defaultdict

import pysam

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    filename2data = dict()
    for filename in args:
        filename2data[filename] = dict(clean_data(get_quast_data(filename)))

    first_time = True
    for key, data in filename2data.items():
        if len(data) < 1:
            continue
        
        if first_time:
            print(",".join(["name", *data.keys()]))
            first_time = False

        print(",".join([key, *data.values()]))
    
        
def get_quast_data(filename):
    quast_data = dict()

    with open(filename) as file_handler:
        reader = csv.reader(file_handler, delimiter='\t')
        for row in reader:
            quast_data[row[0]] = row[1]

    return quast_data

def clean_data(data):
    keeped_key = {"# contigs", "Largest contig", "Total length", "Reference length", "NGA50", "Largest alignment", "# mismatches per 100 kbp", "# indels per 100 kbp", "# misassemblies", "# misassembled contigs"}

    for key, value in data.items():
        if key in keeped_key:
            yield (key, value)

            
if __name__ == "__main__":
    main(sys.argv[1:])
