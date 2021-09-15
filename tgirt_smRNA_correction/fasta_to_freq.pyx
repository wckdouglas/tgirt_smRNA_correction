#!/usr/bin/env python

from itertools import groupby
from collections import defaultdict
from .build_model import make_column_name
import pandas as pd



def fasta_reader(fasta_file):
    """
    fasta iterator:
    adapted from https://www.biostars.org/p/710/
    """


    fasta = open(fasta_file)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    fasta_iterator = (seq_id[1] for seq_id in groupby(fasta, lambda line: line[0] == ">"))
    for seq_record in fasta_iterator:
        # drop the ">"
        seq_id = next(seq_record)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(seq.strip() for seq in next(fasta_iterator))
        yield seq_id, seq


class background_freq:
    '''
    computing the end nucleotide frequencies from a fasta file
    '''
    def __init__(self, fasta_file, num_nucleotide=3):
        self.fasta_file = fasta_file
        self.num_nucleotide = num_nucleotide

        self.end_nucleotide_counter = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.background_nucleotide_df = None
        
        self.compute_end_freq()
        self.compute_table()

    def compute_end_freq(self):
        '''
        count end nucleotide frequencies
        '''
        cdef:
            str seq_id, seq, b5, b3
            int i

        for seq_id, seq in fasta_reader(self.fasta_file):
            for i, (b5, b3) in enumerate(zip(seq[:self.num_nucleotide], seq[-self.num_nucleotide:])):
                self.end_nucleotide_counter["5'"][i][b5] += 1
                self.end_nucleotide_counter["3'"][i][b3] += 1


    def compute_table(self):
        '''
        cmoputing the frequencies and assigning names that matched model
        '''
        ds = []
        for end, end_dict in self.end_nucleotide_counter.items():
            end_label = 'head' if end == "5'" else 'tail'
            ds.append(pd.DataFrame()\
                .from_dict(end_dict, orient='index') \
                .reset_index()
                .pipe(pd.melt, var_name = 'base', value_name = 'base_count', id_vars = 'index') \
                .assign(base_fraction = lambda d: d.groupby(['index'])['base_count']\
                                                    .transform(lambda x: x/x.sum()))\
                .assign(end = end_label))
        self.background_nucleotide_df = pd.concat(ds) \
            .reset_index()\
            .assign(column_name = lambda d: d.end + d['index'].astype(str) + '_' + d.base)\
            .assign(column_name = lambda d: make_column_name(d.column_name, self.num_nucleotide))

    def background(self):
        return self.background_nucleotide_df \
            .filter(['column_name', 'base_fraction'])     \
            .set_index('column_name')

