#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse

import pysam

from collections import defaultdict

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    name2chimera = dict()
    for a in args:
        name2chimera[a] = get_nb_chimera(a)
        
    print("name,nb_chimera")
    for name, cpt in name2chimera.items():
        print(name, cpt)

def get_nb_chimera(filename):
    read2poss = defaultdict(list)
    read2len = defaultdict(int)
    chimera_cpt = 0
    
    mapping = pysam.AlignmentFile(filename, "r")
    for m in mapping:
        if len(m.query_sequence) < 2000:
            continue
        
        if m.infer_query_length() is not None and read2len[m.query_name] < m.infer_query_length():
            read2len[m.query_name] = m.infer_query_length()
            
        read2poss[m.query_name].append((m.reference_name, m.reference_start, m.reference_end))

    for read, map_poss in read2poss.items():
        if len(map_poss) < 2:
            continue
        
        if len(set((m[0] for m in map_poss))) > 1:
            print("read {} is a chimera {}".format(read, map_poss))
            chimera_cpt += 1

        ref_pos = merge_positions(map_poss)
        if read2len[read]/(ref_pos[1] - ref_pos[0]) < 0.9:
            print("read {} is a chimera {}".format(read, read2len[read]/(ref_pos[1] - ref_pos[0])))
            chimera_cpt += 1
            
    return chimera_cpt

def merge_positions(map_poss):
    base = list(map_poss.pop(0)[1:])

    for (_, beg, end) in map_poss:
        if beg < base[0]:
            base[0] = beg
        if end > base[1]:
            base[1] = end
        
    return base

import math
def get_sequence_entropy(seq):
    e = 0
    for n in "ACTG":
        x = seq.count(n)/len(seq)
        if x > 0:
            e -= x * math.log(x, 2)

    return e

        
if __name__ == "__main__":
    main(sys.argv[1:])
