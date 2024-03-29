from __future__ import print_function
from itertools import groupby
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
from libc.math cimport ceil, log2, exp
import csv
from xopen import xopen
from .model import h2o_rf
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__file__)

def total_fragment_count(bed_file):
    '''
    Count total number of fragments
    '''
    with open(bed_file,'r') as inbed:
        total = sum(1 for line in inbed)
    return total

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
    seq = fa.fetch(str(chrom), int(start), int(end))
    seq = seq if strand == "+" else reverse_complement(seq)
    return seq.upper()


def seq_to_bias(bias_idx, fa, chrom,start,end, strand):
    '''
    feed in coordinates, fetch sequence, then fetch bias factor
    '''

    seq = fetch_seq(fa, chrom, start, end, strand)

    return extract_correction_factor(bias_idx,  seq)



def feature_engineering(d):
    ds = []
    for end in ['head', 'tail']:
        ds.append(pd.DataFrame(d['%s_seq' %end]\
                    .apply(list)\
                    .tolist(), 
                columns = ['%s%i' %(end,i) for i in range(3)]))
    features = pd.concat(ds, axis=1)
    return pd.concat([d, features], axis=1) 
    

def check_cols(df, cols):
    existing_cols = set(df.columns)
    required_cols = set(cols)
    need_add_cols = required_cols - existing_cols
    for col in need_add_cols:
        df[col] = 0
    return df


def model_correction(bed, genome_fa, bias_index, outfile):
    """
    :param str bed: bed file to be weigted
    :param pysam.FastaFile genome_fa: genome fasta file that the bed file mapped to
    :param bias_index: dictionary of weights or model
    :param filehandle outfile: outfile path
    """
    # genome_fa = pysam.Fastafile'/stor/work/Lambowitz/ref/RNASeqConsortium/ercc/ERCC92.fa')
    if 'model' in bias_index:
        use_rf_model(bed, genome_fa, bias_index, outfile)
    else:
        use_weight_model(bed, genome_fa, bias_index, outfile)

def use_weight_model(bed, genome_fa, bias_index, outfile):
    cdef:
        str line

    logger.info('Using random forest model')
    with xopen(bed) as inbed:
        for line in inbed:
            chrom, start, end, strand = itemgetter(0,1,2,5)(line.split('\t'))
            seq = fetch_seq(genome_fa, chrom, start, end, strand)
            weight = bias_index[seq[:3] + ',' + reverse_complement(seq[-3:])]
            print('{}\t{}'.format(line.strip(),weight), file=outfile)



def use_rf_model(bed, genome_fa, bias_index, outfile):
    logger.info('Using random forest model')
    model_path = bias_index['model']
    cols = bias_index['X_col']
    model = h2o_rf()
    model.load_model(model_path)
    get_seq = partial(fetch_seq, genome_fa)

    for i, bed_df in enumerate(pd.read_table(bed,header=None, chunksize=1000000)):
        out_cols = bed_df.columns.tolist()
        out_cols = list(map(lambda x: 'X%i' %x, out_cols))
        bed_df.columns = out_cols
        out_cols.append('correction_factor')
        
        out_df = bed_df \
            .assign(seq = lambda d: list(map(get_seq, d['X0'], d['X1'], d['X2'], d['X5'])))\
            .assign(head_seq = lambda d: d.seq.str.slice(0,3))\
            .assign(tail_seq = lambda d: (d.seq + 'N').str.slice(-4,-1)) \
            .reset_index() \
            .pipe(feature_engineering)  \
            .pipe(check_cols, cols)\
            .assign(pred = lambda d: model.predict(d.loc[:, cols]))  \
            .assign(nucleotide = lambda d: d.loc[:, cols].sum(axis=1))\
            .assign(correction_factor = lambda d: np.exp(-d.pred))  \
            .assign(correction_factor = lambda d: d.shape[0] * d.correction_factor/d.correction_factor.sum())

        #assert(out_df.query('nucleotide != 6').shape[0] == 0)
        
        mode = 'w' if i ==0 else 'a'
        out_df\
            .loc[:, out_cols]\
            .to_csv(outfile, mode = mode, header=False, sep='\t', index=False, float_format='%.6f')


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
        for coor, lines in groupby(inbed, group_func):
            chrom, start, end, strand = coor
            lines = list(lines)
            bf = fetch_bf_func(chrom, int(start), int(end), strand)
            bf = len(lines) * bf
            assert(bf > 0)
            print('{}\t{}'.format(lines[0].strip(), bf), file=outfile)
       






        

