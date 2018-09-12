
from __future__ import print_function
from itertools import groupby
from libc.math cimport ceil, log2, exp
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

def fetch_seq(fa, chrom, start, end, strand):
    seq = fa.fetch(chrom, start, end)
    seq = seq if strand == "+" else reverse_complement(seq)
    return seq


def seq_to_bias(bias_idx, fa, chrom,start,end, strand):
    '''
    feed in coordinates, fetch sequence, then fetch bias factor
    '''

    seq = fetch_seq(fa, chrom, start, end, strand)

    return extract_correction_factor(bias_idx,  seq)


def parse_bed(bed, genome_fa, bias_index, outfile):
    '''
    for each bed line, extract sequence and first 3 nucleotides from each end
    and extract bias factor
    '''
    cdef:
        str line, chrom, strand, name
        str start, end
        int number_of_out
        int out_count = 0
        int in_count = 0
        double bf, log_out_count
        int highest_cov = 0,
        int line_count
        int i
        double pred_count

    fetch_bf_func = partial(seq_to_bias, bias_index, genome_fa)

    with open(bed) as inbed:
        for line in inbed:
            fields = line.strip().split('\t')
            chrom, start, end, strand = itemgetter(0,1,2,5)(fields)
            bf = fetch_bf_func(chrom, int(start), int(end), strand)
            assert(bf > 0)
            print('{}\t{}'.format(line.strip(), bf), file=outfile)
       






        


