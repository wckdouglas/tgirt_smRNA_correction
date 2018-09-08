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


def seq_to_bias(bias_idx, fa, chrom,start,end, strand):
    '''
    feed in coordinates, fetch sequence, then fetch bias factor
    '''
    seq = fa.fetch(chrom, start, end)
    seq = seq if strand == "+" else reverse_complement(seq)

    return extract_correction_factor(bias_idx,  seq)


def parse_bed(bed, genome_fa, index_file, outfile):
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

    with open(index_file,'rb') as idx:
        bias_index = pickle.load(idx) 
    

    fetch_bf_func = partial(seq_to_bias, bias_index, genome_fa)
    bed_template = '{chrom}\t{start}\t{end}\tFrag_{frag_count}_{line_count}\t{frag_len}\t{strand}'
    rows = []
    with open(bed) as inbed:
        for coor, bedlines in groupby(inbed, key=group_func):
            chrom, start, end, strand = coor
            line_count = len(list(bedlines))
            rows.append((chrom, start, end, strand, line_count))
            in_count += line_count
        
    print ('Read %i fragments' %in_count)
    
    bed_df = pd.DataFrame(rows, names = ['chrom','start','end','strand','RC']) \
        .assign(RC_frac = lambda d: d.RC/in_count) \
        .assign(bf = lambda d: list(map(fetch_bf_func, d.chrom, d.start.astype(int),
                                        d.end.astype(int), d.strand))) \
        .assign(log2_fraction = lambda d: np.log2(RC_frac) - d.bf ) \
        .assign(out_fraction = lambda d: d.log2_fraction.rpow(2)) \
        .assign(out_count = lambda d: d.out_fraction * in_count) \
        .query('out_count >= 1') \
        .assign(out_count = lambda d: d.out_count.astype(int)) \
        .assign(frag_length = lambda d: d.end.astype(int) - d.start.astype(int) )
        
    for line_count, row in bed_df.iterrows():
        for i in range(row['out_count']):
            out_line = bed_template.format(chrom = row['chrom'],
                                            start = row['start'],
                                            end = row['end'],
                                            frag_count = line_count,
                                            line_count = i,
                                            frag_len = row['frag_length'],
                                            strand = row['strand]'')
            print(out_line, file = outfile)
            out_count += 1

    print('Parsed: ', in_count, '\nPrinted: ', out_count, file = sys.stderr)






        



