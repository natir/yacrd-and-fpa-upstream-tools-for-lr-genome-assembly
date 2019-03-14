#!/usr/bin/env python3

## INSPIRED BY DASCRUBBER

import sys
import argparse
import itertools

def rename_reads_with_fake_pacbio_names(read_filein, read_fileout):

    read_name_dict, read_comment_dict = {}, {}
    read_names = set()

    counter = itertools.count(start=0)
    base_total = 0

    with open(read_filein, 'rt') as seq_file:
        with open(read_fileout, 'wt') as renamed_reads:
            read_num, update_interval = 0, 1

            for header in seq_file:
                try:
                    assert header[0] == '>'
                    header = header[1:].strip()
                    assert len(header) > 0
                except AssertionError:
                    sys.exit('Error: failed to parse read header')

                header_parts = header.split(' ', 1)
                old_name = header_parts[0]
                if old_name in read_names:
                    sys.exit('Error: duplicate read name: ' + old_name)
                read_names.add(old_name)
                try:
                    old_comment = header_parts[1]
                except IndexError:
                    old_comment = None

                seq = next(seq_file).strip()
                read_num = next(counter)
                new_header = 'reads/' + str(read_num) + '/0_' + str(len(seq))

                read_name_dict[read_num] = old_name
                read_comment_dict[read_num] = old_comment

                renamed_reads.write('>')
                renamed_reads.write(new_header)
                renamed_reads.write('\n')
                renamed_reads.write(seq)
                renamed_reads.write('\n')

                base_total += len(seq)

def int_to_str(num):
    return '{:,}'.format(num)

if __name__ == "__main__":
    rename_reads_with_fake_pacbio_names(sys.argv[1], sys.argv[2])
