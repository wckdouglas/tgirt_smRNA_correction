#!/usr/bin/env python

from __future__ import print_function
import pysam
import argparse
import pyximport
from tgirt_smRNA_correction.bed_parser import parse_bed
import sys

def get_opt():
    parser = argparse.ArgumentParser(description='Using prebuilt bias index, add a column to bed fragment file as correction factor')
    parser.add_argument('-f','--fasta', help='Genome fasta file', required=True)
    parser.add_argument('-b','--bed', help='Input fragment bed file (default: -)', default='-')
    parser.add_argument('-i','--index', help = 'Index', required=True)
    parser.add_argument('-o','--out_bed', default='-',help='output file (default: -)')
    args = parser.parse_args()
    return args

def main():
    args = get_opt()
    fa = pysam.FastaFile(args.fasta)
    bed = sys.stdin if args.bed == "-" or args.bed == '/dev/stdin' else args.bed
    outfile = sys.stdout if args.out_bed == '-' or args.out_bed == '/dev/stdout' else open(args.out_bed,'w')
    parse_bed(bed, fa, args.index, outfile)


if __name__ == '__main__':
    main()
