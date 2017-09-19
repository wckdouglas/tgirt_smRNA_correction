#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np
import string
from collections import Counter
import argparse
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder, StandardScaler


def fasta_parser(fasta_file_handle):
    seq = ''
    seq_id = ''
    while True:
        line = fasta_file_handle.next()
        if line.startswith('>'):
            if seq_id and seq:
                yield seq_id, seq
            seq_id = line.lstrip('>').rstrip('\n').split(' ')[0]
            seq = ''
        else:
            seq += line.strip('\n')
    yield seq_id, seq


def make_header(nucleotide):
    head_title = '\t'.join(map(lambda x: 'head%i' %x, range(nucleotide)))
    tail_title = '\t'.join(map(lambda x: 'tail%i' %x, range(nucleotide)))
    return 'seq_id\t%s\t%s' %(head_title,tail_title)


def GC_content(seq):
    base_content = Counter(seq)
    return (base_content['G'] + base_content['C'])/float(len(seq))


class input_table():

    def __init__(self, args):
        self.base_table = args.output_prefix + '.base_table.tsv'
        self.fasta_file = args.fasta
        self.nucleotide_count = args.nucleotide
        self.expected_count_file = args.expected_count
        self.experimental_count_file = args.experimental_count
        self.out_table_name = args.output_prefix + '_train_set.tsv'
        self.out_table = None


    def build_base_table(self):
        with open(self.fasta_file,'r') as fasta, open(self.base_table,'w') as out_file :
            print(make_header(self.nucleotide_count), file=out_file)
            for seq_id, seq in fasta_parser(fasta):
                head_3_base = '\t'.join(list(seq[:self.nucleotide_count]))
                tail_3_base = '\t'.join(list(seq[-self.nucleotide_count:]))
                print('%s\t%s\t%s' %(seq_id, head_3_base, tail_3_base), file=out_file)

        self.out_table =  pd.read_table(self.base_table)

    def incorporate_count(self, input_table, column_name):
        if input_table:
            df = pd.read_csv(input_table,
                         names = ['seq_id',column_name]) \
                .merge(self.out_table, how='right')\
                .fillna(0)

            if df.shape[0] < self.out_table.shape[0]:
                print( '%i sequence not in %s' %(self.out_table.shape[0]-df.shape[0], column_name), file=sys.stderr)
        else:
            self.out_table[column_name] = 10000
            df = self.out_table
        
        self.out_table = df

    def add_observations(self):
        self.incorporate_count(self.expected_count_file,'expected_count')
        self.incorporate_count(self.experimental_count_file ,'experimental_count')

    def make_table(self):
        self.out_table.to_csv(self.out_table_name, sep='\t', index=False)
        print('Written:', self.out_table_name, file = sys.stderr)

