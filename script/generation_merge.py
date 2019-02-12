#!/usr/bin/env python3

import sys
import gzip
import random
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(args):

    if args == None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, help="read input",
                        type=argparse.FileType('r'))
    parser.add_argument("-o", "--output", required=True, help="read output",
                        type=argparse.FileType('w'))
    parser.add_argument("-n", "--number", default=1000,
                        help="number of read generate", type=int)
    parser.add_argument("-s", "--seed", default=42,
                        help="set the random seed default 42", type=int)

    args = vars(parser.parse_args(args))
    
    random.seed(args["seed"])

    reader = SeqIO.parse(args["input"], "fasta")

    for i, (read1, read2) in enumerate(zip(reader, reader)):
        if args["number"] > 0 and i >= args["number"]:
            write_read(read1.id, read1.seq, args["output"])
            write_read(read2.id, read2.seq, args["output"])
            break

        insert = "".join(random.sample(['A', 'C', 'T', 'G']*13, 50))
        write_read("chimeric_"+read1.id+"_"+read2.id,
                   str(read1.seq) + insert + str(read2.seq),
                   args["output"])

    for i in reader:
        write_read(i.id, i.seq, args["output"])


def write_read(iden, seq, output):
    print(">" + iden, file=output)
    print(seq, file=output)

if __name__ == "__main__":
    main(sys.argv[1:])
