#!/usr/bin/env python
# coding=utf-8

import argparse
import os
import sys

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def consensus(aln):
    # precompute dictionaries on first run
    if "collapse_iupac" not in consensus.__dict__:
        # sorted tuples map to iupac codes
        consensus.collapse_iupac = {
            ('-',): '-',
            ('a',): 'a',
            ('g',): 'g',
            ('c',): 'c',
            ('t',): 't',
            ('c', 't'): 'y',
            ('a', 'g'): 'r',
            ('a', 't'): 'w',
            ('c', 'g'): 's',
            ('g', 't'): 'k',
            ('a', 'c'): 'm',
            ('a', 'g', 't'): 'd',
            ('a', 'c', 'g'): 'v',
            ('a', 'c', 't'): 'h',
            ('c', 'g', 't'): 'b',
            ('a', 'c', 'g', 't'): 'n',
        }

        # inverse process of collapse_iupac, but n's in input are treated special, so don't expand them
        consensus.expand_iupac = {value: set(key) for key, value in consensus.collapse_iupac.items() if value != 'n'}
        consensus.expand_iupac['n'] = {'n'}

    con = ""
    for loc in range(aln.get_alignment_length()):
        # expand iupac code to base set
        base_set = set.union(*[consensus.expand_iupac[seq[loc].lower()] for seq in aln])

        # if vertical '-', return '-'
        if base_set == {'-'}:
            con += '-'
            continue

        # if vertical 'n', return 'n'
        if base_set == {'n'}:
            con += 'n'
            continue

        # if not vertical '-', ignore '-'
        if '-' in base_set:
            base_set.remove('-')

        # if not vertical 'n', ignore 'n'
        if 'n' in base_set:
            base_set.remove('n')

        # collapse base_set to iupac code
        con += consensus.collapse_iupac[tuple(sorted(base_set))]

    return con


class HelpAndQuitOnFailParser(argparse.ArgumentParser):
    """custom argparse configuration
    if error parsing, prints help and exits"""

    def error(self, message):
        sys.stderr.write('error: {}\n'.format(message))
        self.print_help()
        sys.exit(2)


def main():
    parser = HelpAndQuitOnFailParser()

    # files/directories
    parser.add_argument('-i', '--in_directory', default="input",
                        help='path to input directory')

    parser.add_argument('-o', '--out_file', default="consensus.fasta",
                        help='path to output file')

    args = parser.parse_args()

    with open(args.out_file, "wt") as f:
        for filename in os.listdir(args.in_directory):

            sample_name, extension = os.path.splitext(filename)
            if extension not in {".fasta", ".fna"}:
                continue

            in_aln = AlignIO.read(os.path.join(args.in_directory, filename), "fasta")

            con_str = consensus(in_aln)
            con_seq = Seq(con_str, alphabet=IUPAC.ambiguous_dna)
            con_rec = SeqRecord(con_seq, id=sample_name, name=sample_name, description="")

            f.write(con_rec.format("fasta"))


if __name__ == "__main__":
    # run main
    sys.exit(main())
