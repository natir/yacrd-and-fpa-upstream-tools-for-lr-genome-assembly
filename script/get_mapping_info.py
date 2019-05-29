#!/usr/bin/env python3

import os
import csv
import sys

from collections import defaultdict

import pysam


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    filename2data = dict()
    for filename in args:
        nb_map, edit_distance_sum, mapping_length, s_lengths = parse_mapping(filename)
        n10, l10 = get_N_L(s_lengths, 0.1)
        n50, l50 = get_N_L(s_lengths, 0.5)
        n90, l90 = get_N_L(s_lengths, 0.9)

        filename2data[filename] = {
            "nb read map": str(nb_map),
            "edit distance": str(edit_distance_sum),
            "number of read": str(len(s_lengths)),
            "mapping len": str(mapping_length),
            "n10": str(n10),
            "n50": str(n50),
            "n90": str(n90),
            "l10": str(l10),
            "l50": str(l50),
            "l90": str(l90),
        }

        
    first_time = True
    for key, data in filename2data.items():
        if len(data) < 1:
            continue
        
        if first_time:
            print(",".join(["name", *data.keys()]))
            first_time = False

        print(",".join([key, *data.values()]))
    

def parse_mapping(mapping_file):
    read_id = set()
    read_len = list()
    edit_distance = 0
    mapping_length = 0
    
    if not os.path.isfile(mapping_file):
        return -1, -1
    
    mapping = pysam.AlignmentFile(mapping_file, "rb")
    for m in mapping:
        if m.query_length < 2000:
            continue
        
        if m.flag == 0 or m.flag == 16:
            read_id.add(m.query_name)
            read_len.append(m.infer_query_length())
            edit_distance += m.get_tag("NM")
            mapping_length += sum([count for (t, count) in m.cigartuples if t == 0])
            
    return str(len(read_id)), str(edit_distance), str(mapping_length), read_len


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
    
    return str(index), str(lengths[index])


if __name__ == "__main__":
    main(sys.argv[1:])
