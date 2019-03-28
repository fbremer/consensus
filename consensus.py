#!/usr/bin/env python
# coding=utf-8

import argparse
import os
import sys

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def consensus(aln, any_n=False, any_gap=False):
    # pre-compute dictionaries on first run
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

        # with any_gap, any '-' in input will result in '-' in output
        if any_gap and '-' in base_set:
            con += '-'
            continue
        # without any_gap, if vertical '-', return '-'
        elif not any_gap and base_set == {'-'}:
            con += '-'
            continue
        # if not vertical '-', ignore '-'
        elif '-' in base_set:
            base_set.remove('-')

        # with any_n, any 'n' in input will result in 'n' in output
        if any_n and 'n' in base_set:
            con += 'n'
            continue
        # without any_n, if vertical 'n', return 'n'
        elif not any_n and base_set == {'n'}:
            con += 'n'
            continue
        # if not vertical 'n', ignore 'n'
        elif 'n' in base_set:
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
    parser.add_argument('-i', '--input', default="alignment.fasta",
                        help='path to input file, or input directory if using batch mode')

    parser.add_argument('-o', '--out_file', default="consensus.fasta",
                        help='path to output file')

    parser.add_argument('--batch_mode', action='store_true',
                        help='take an input directory instead of input file')

    parser.add_argument('--any_n', action='store_true',
                        help=("By default, if *all* bases at a location are 'n', the consensus will be a 'n'. "
                              "This flag changes that rule to return an 'n' if *any* 'n' are found at that loc."))

    parser.add_argument('--any_gap', action='store_true',
                        help=("Same as --any_n above. "
                              "If --any_gap and --any_n are set, gaps take precedence"))

    args = parser.parse_args()

    with open(args.out_file, "wt") as f:
        # batch mode
        if args.batch_mode:
            for filename in os.listdir(args.input):

                sample_name, extension = os.path.splitext(filename)
                if extension not in {".fasta", ".fna"}:
                    continue

                in_aln = AlignIO.read(os.path.join(args.input, filename), "fasta")

                con_str = consensus(in_aln, any_n=args.any_n, any_gap=args.any_gap)
                con_seq = Seq(con_str, alphabet=IUPAC.ambiguous_dna)
                con_rec = SeqRecord(con_seq, id=sample_name, name=sample_name, description="")

                f.write(con_rec.format("fasta"))

        # single file mode
        else:
            sample_name, _ = os.path.splitext(args.input)

            in_aln = AlignIO.read(args.input, "fasta")

            con_str = consensus(in_aln, any_n=args.any_n, any_gap=args.any_gap)
            con_seq = Seq(con_str, alphabet=IUPAC.ambiguous_dna)
            con_rec = SeqRecord(con_seq, id=sample_name, name=sample_name, description="")

            f.write(con_rec.format("fasta"))


if __name__ == "__main__":
    # run main
    sys.exit(main())
