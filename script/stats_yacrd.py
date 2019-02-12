#!/usr/bin/env python3

import sys
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(args):

    if args == None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--chimeric", required=True, help="chimeric read",
                        type=argparse.FileType('r'))
    parser.add_argument("-y", "--yacrd", required=True, help="yacrd output",
                        type=argparse.FileType('r'))


    args = vars(parser.parse_args(args))

    true_chimeric = set()
    true_not_chimeric = set()
    for read in SeqIO.parse(args["chimeric"], "fasta"):
        if "chimeric" in read.id:
            true_chimeric.add(read.id)
        else:
            true_not_chimeric.add(read.id)

    all_read = true_chimeric | true_not_chimeric

    positif_notcov = set()
    positif_chimeric = set()
    for row in csv.reader(args["yacrd"], delimiter="\t"):
        if row[0] == "Not_covered":
            positif_notcov.add(row[1])
        if row[0] == "Chimeric":
            positif_chimeric.add(row[1])

    negatif = all_read - positif_chimeric
    
    P = len(true_chimeric)
    N = len(true_not_chimeric)

    FN = len(true_chimeric & negatif)
    TN = len(true_not_chimeric & negatif)
    TP = len(true_chimeric & positif_chimeric)
    FP = len(true_not_chimeric & positif_chimeric)

    print("positif", len(true_chimeric))
    print("negatif", len(true_not_chimeric))
    try:
        print("precision : {:.2%}".format(TP / (TP + FP)))
    except ZeroDivisionError:
        print("can't compute precision TP: {} (TP + FP): {}".format(TP, TP + FP))
    
    try:
        print("sensitivity : {:.2%}".format(TP / P))
    except ZeroDivisionError:
        print("can't compute sensitivity TP: {} P: {}".format(TP, P))
    
    try:
        print("specificity : {:.2%}".format(TN / N))
    except ZeroDivisionError:
        print("can't compute specificity TN: {} N: {}".format(TN, N))

if __name__ == "__main__":
    main(sys.argv[1:])
