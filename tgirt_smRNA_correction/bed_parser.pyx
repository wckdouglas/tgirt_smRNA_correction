from __future__ import print_function
from itertools import groupby
from libc.math cimport ceil
from scipy.stats import poisson
import pysam
from operator import itemgetter
import string
import sys
import pickle
import six
from builtins import range, map
from functools import partial
import pandas as pd
import numpy as np

def total_fragment_count(bed_file):
    '''
    Count total number of fragments
    '''
    with open(bed_file,'r') as inbed:
        total = sum(1 for line in inbed)
    return total

try:
    complement_seq = string.maketrans('ACTGNactgn','TGACNtgacn')
except AttributeError:
    complement_seq = str.maketrans('ACTGNactgn','TGACNtgacn')

def complement(seq):
    """
    Find complement a sequence.
    ============================
    parameter:
    
    string seq: sequence to be complemented
    
    return:
    complemented sequence
    """
    return seq.translate(complement_seq)


def reverse_complement(seq):
    """
    Reverse complement a sequence.
    ============================
    parameter:
    
    string seq: sequence to be reverse complemented
    
    return:
    reverse complemented sequence
    """
    return complement(seq)[::-1]

def group_func(x):
    fields = x.rstrip().split('\t')
    return itemgetter(0,1,2,5)(fields)


def extract_correction_factor(bias_idx, seq):
    '''
    given a sequence, extract head and tail and pull bias factor
    '''
    key = seq[:3] + ','+ seq[-3:]
    correction_factor = bias_idx[key.upper()]
    return correction_factor


def seq_to_bias(bias_idx, fa, chrom,start,end, strand):
    '''
    feed in coordinates, fetch sequence, then fetch bias factor
    '''
    seq = fa.fetch(chrom, start, end)
    seq = seq if strand == "+" else reverse_complement(seq)

    return extract_correction_factor(bias_idx, seq)


def make_bed_df(bed):
    '''
    input a bed file, groupby chrom,start,end ,strand and add a count value to these coordinat
    '''
    rows = []
    with open(bed, 'r') as inbed:
        for coor, bedlines in groupby(inbed, group_func):
            count = len(list(bedlines))
            coor = list(coor)
            coor.append(count)
            rows.append(coor)

    return pd.DataFrame(rows, columns = ['chrom','start','end', 'strand', 'count']) \
        .assign(cpm = lambda d: d['count'].transform(lambda x: x/x.sum()))  \
        .assign(log2cpm = lambda d: np.log2(d['cpm']))


def output_bed(bed_df, outfile):
    '''
    given a corrected count bed, output bed lines
    '''

    template = '{chrom}\t{start}\t{end}\tFrag_{frag_count}\t{frag_len}\t{strand}'

    for i, (idx, bedline) in enumerate(bed_df.iterrows()):
        line = template.format(chrom = bedline['chrom'],
                                start = bedline['start'],
                                end = bedline['end'],
                                frag_count = i, 
                                strand = bedline['strand'])
        print(line, file = outfile)


def parse_bed(bed, genome_fa, index_file, outfile):
    '''
    for each bed line, extract sequence and first 3 nucleotides from each end
    and extract bias factor
    '''
    cdef:
        str line, chrom, strand, name
        int start, end
        int number_of_out
        int out_count = 0

    with open(index_file,'rb') as idx:
        bias_index = pickle.load(idx) 

    fetch_bf_func = partial(seq_to_bias, bias_index, genome_fa)
    bed_df = make_bed_df(bed) \
        .assign(correction_factor = lambda d: list(map(fetch_bf_func, d.chrom, d.start.astype(int), d.end.astype(int), d.strand))) \
        .assign(corrected_log2_cpm = lambda d: d.log2cpm - d.correction_factor)\
        .assign(pseudo_count = lambda d: d.corrected_log2_cpm.rpow(2))  \
        .assign(new_count = lambda d: d['count'].min() * d.pseudo_count/d.pseudo_count.min())

    
    output_bed(bed_df, outfile)




    





        



