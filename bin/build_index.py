#!/usr/bin/env python

import argparse
from tgirt_smRNA_correction.build_table import input_table
from tgirt_smRNA_correction.build_model import lm_model

def get_opt():
    parser = argparse.ArgumentParser(description='Building nucleotide bias table for each transcripts')
    parser.add_argument('-f','--fasta', help='Small RNA fasta file', required=True)
    parser.add_argument('-n','--nucleotide', help='How many nucleotide to look at from both end?', default=3, type=int)
    parser.add_argument('-o','--output_prefix', default='-',help='How many nucleotide to look at from both end?')
    parser.add_argument('-e','--expected_count', help='Expected count (comma delimintaed: seq,count) (default: all equal)')
    parser.add_argument('-c','--experimental_count', help='Experimental count (comma delimintaed: seq,count)', required=True)
    args = parser.parse_args()
    return args

def main():
    args = get_opt()
    
    ## make nucleotide table
    base_table = input_table(args)
    base_table.build_base_table()
    base_table.add_observations()
    base_table.make_table()
    train_table = base_table.out_table_name


    ## train model
    index_file = args.output_prefix + '_index.pkl'
    lm = lm_model(train_table, index_file)
    lm.preprocess_data()
    lm.train_lm()
    lm.write_index()

if __name__ == '__main__':
    main()
