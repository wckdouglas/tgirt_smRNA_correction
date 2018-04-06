from __future__ import print_function
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

def parse_bed_line(str bed_line):
    '''
    output chrom, start, end, strand, name fields from bed line
    '''
    cdef:
        str chrom, start, end, strand, name
    
    fields = bed_line.strip().split('\t')
    chrom, start, end, name, strand = itemgetter(0,1,2,3,5)(fields)
    return chrom, int(start), int(end), name, strand

cdef int extract_correct_factor(dict bias_index, float observed_cpm, str seq):
    '''
    correction factor is cpm offset

    original_cpm =  1/total * 1e6
    log2_original_cpm = log2(original_cpm)
    corrected_log2_cpm = log2_original_cpm - log2(correction_factor)
    corrected_cpm = 2^corrected_cpm


    ---
    can be reduced to:
        corrected = obs_cpm/(correction_factor)
    '''
    cdef:
        str key
        float correction_factor
        double corrected_cpm
        int corrected_count

    key = seq[:3] + ','+ seq[-3:]
    correction_factor = bias_index[key.upper()]
    corrected_cpm = ceil(observed_cpm/(correction_factor))
    corrected_count = poisson.rvs(corrected_cpm)

    return  corrected_count


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


    total_fragment = total_fragment_count(bed)
    observed_cpm = 1e6 / total_fragment 
    extract_correction = partial(extract_correct_factor, bias_index, observed_cpm)

    for line in open(bed, 'r'):
        chrom, start, end, name, strand = parse_bed_line(line)
        seq = genome_fa.fetch(chrom, start, end)
        if 'N' not in seq:
            seq = reverse_complement(seq) if strand == "-" else seq
            number_of_out = extract_correction(seq)
            for i in range(number_of_out):
                print(line.strip(), file=outfile)
                out_count += 1
        else:
            print('Skipping: %s sequence contain base N' %(name),file=sys.stderr)
    
    print('Output %i fragments ' %(out_count), file=sys.stderr)
