#!/usr/bin/env python

from __future__ import print_function
import pysam
import argparse
import pyximport
from tgirt_smRNA_correction.bed_parser import parse_bed

def get_opt():
    parser = argparse.ArgumentParser(description='Building nucleotide table for each transcripts')
    parser.add_argument('-f','--fasta', help='Genome fasta file', required=True)
    parser.add_argument('-b','--bed', help='Input fragment bed file (default: -)', default='-')
    parser.add_argument('-i','--index', help = 'Index', default='../test/index.pkl')
    parser.add_argument('-o','--output_prefix', default='-',help='')
    args = parser.parse_args()
    return args

def main():
    args = get_opt()
    fa = pysam.FastaFile(args.fasta)
    parse_bed(args.bed, fa, args.index)


if __name__ == '__main__':
    main()
