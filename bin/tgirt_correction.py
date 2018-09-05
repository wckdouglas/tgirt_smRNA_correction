#!/usr/bin/env python

from __future__ import print_function
import pysam
import argparse
import pyximport
from tgirt_smRNA_correction.bed_parser import parse_bed
from tgirt_smRNA_correction.build_table import input_table
from tgirt_smRNA_correction.build_model import lm_model
import sys
import os

def get_opt():
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]))
    subparser = parser.add_subparsers(help = 'Run type:', dest = 'subcommand')
    subparser.required=True
        
    #add build index
    index_builder = subparser.add_parser('build', 
                    help='Building nucleotide bias table for each transcripts')
    index_builder.add_argument('-f','--fasta', help='Small RNA fasta file', required=True)
    index_builder.add_argument('-n','--nucleotide', help='How many nucleotide to look at from both end?', default=3, type=int)
    index_builder.add_argument('-o','--output_prefix', default='-',help='How many nucleotide to look at from both end?')
    index_builder.add_argument('-e','--expected_count', help='Expected count (comma delimintaed: seq,count) (default: all equal)')
    index_builder.add_argument('-c','--experimental_count', help='Experimental count (comma delimintaed: seq,count)', required=True)


    
    correction = subparser.add_parser('correct', 
            help='Using prebuilt bias index, add a column to bed fragment file '\
                'as log2(CPM_observed) - log2(CPM_expected)')
    
    correction.add_argument('-f','--fasta', help='Genome fasta file', required=True)
    correction.add_argument('-b','--bed', help='Input fragment bed file', required=True)
    correction.add_argument('-i','--index', help = 'Index', required=True)
    correction.add_argument('-o','--out_bed', default='-',help='output file (default: -)')

    args = parser.parse_args()
    return args


def build(args):
    # building index    
    ## make nucleotide table
    base_table = input_table(args)
    base_table.build_base_table()
    base_table.add_observations()
    base_table.make_table()
    train_table = base_table.out_table_name


    ## train model
    index_file = args.output_prefix + '_index.pkl'
    model_file = args.output_prefix + '_model.pkl'
    lm = lm_model(train_table, index_file)
    lm.preprocess_data()
    lm.train_lm()
    lm.write_index()


def correction(args):
    fa = pysam.FastaFile(args.fasta)
    outfile = sys.stdout if args.out_bed == '-' or args.out_bed == '/dev/stdout' else open(args.out_bed,'w')
    parse_bed(args.bed, fa, args.index, outfile)


def main():
    args = get_opt()

    if args.subcommand == 'build':
        build(args)
    
    elif args.subcommand == 'correct':
        correction(args)



if __name__ == '__main__':
    main()
