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

    res_df = pd.DataFrame(rows, columns = ['chrom','start','end', 'strand', 'count']) \
        .groupby(['chrom','start','end','strand'], as_index=False)\
        .agg({'count':'sum'})\
        .assign(cpm = lambda d: d['count'].transform(lambda x: x/x.sum()))  \
        .assign(log2cpm = lambda d: np.log2(d['cpm'])) \
        .assign(start = lambda d: d.start.astype(int)) \
        .assign(end = lambda d: d.end.astype(int))
    print('Read in %i fragments ' %res_df['count'].sum(), file=sys.stderr)
    return res_df


def output_bed(bed_df, outfile):
    '''
    given a corrected count bed, output bed lines
    '''

    template = '{chrom}\t{start}\t{end}\tFrag_{frag_count}_{line_count}\t{frag_len}\t{strand}'
    out_count = 0

    for i, (idx, bedline) in enumerate(bed_df.iterrows()):
        for line_count in range(int(bedline['new_count'])):
            line = template.format(chrom = bedline['chrom'],
                                    start = bedline['start'],
                                    end = bedline['end'],
                                    frag_count = i, 
                                    line_count = line_count, 
                                    frag_len = bedline['end'] - bedline['start'],
                                    strand = bedline['strand'])
            print(line, file = outfile)
            out_count += 1

    print('Output in %i fragments' %(out_count))


def parse_bed_df(bed, genome_fa, index_file, outfile):
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
        .assign(correction_factor = lambda d: list(map(fetch_bf_func, d.chrom, d.start, d.end, d.strand))) \
        .assign(corrected_log2_cpm = lambda d: d.log2cpm - d.correction_factor)\
        .assign(pseudo_count = lambda d: d.corrected_log2_cpm.rpow(2))  \
        .assign(new_count = lambda d: d['count'].sum() * d.pseudo_count/d.pseudo_count.sum()) \
        .assign(new_count = lambda d: d['new_count'].astype(int))  \
        .query('new_count > 0')
    
    bed_df.to_csv('/stor/home/cdw2854/miRNA_seq/new_NTT_corrections/temp.csv', index=False)
    output_bed(bed_df, outfile)



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
    with open(bed) as inbed:
        for coor, bedlines in groupby(inbed, key=group_func):
            chrom, start, end, strand = coor
            line_count = len(list(bedlines))
            bf = fetch_bf_func(chrom, int(start), int(end), strand)
            log2_pred_count = log2(line_count) - bf
            pred_count = 2**(log2_pred_count)

            for i in range(int(pred_count)):
                out_line = bed_template.format(chrom = chrom,
                                                start = start,
                                                end = end,
                                                frag_count = line_count,
                                                line_count = i,
                                                frag_len = int(end)- int(start),
                                                strand = strand)
                print(out_line, file = outfile)

            # update variables
            highest_cov = max(highest_cov, line_count)
            in_count += line_count
            out_count += int(pred_count)
    
    print('Parsed: ', in_count, '\nPrinted: ', out_count, file = sys.stderr)
    if highest_cov == 1:
        sys.exit('Is the input bed sorted?')






        



