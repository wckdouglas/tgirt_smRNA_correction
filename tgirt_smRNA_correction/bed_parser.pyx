#!/bin/env python

from __future__ import print_function
import pysam
from operator import itemgetter
import string
import sys
import pickle


complement_seq = string.maketrans('ACTGNactgn','TGACNtgacn')
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
    cdef:
        str chrom, start, end, strand
    
    fields = bed_line.strip().split('\t')
    chrom, start, end, strand = itemgetter(0,1,2,5)(fields)
    return chrom, int(start), int(end), strand

def parse_bed(bed_file, genome_fa, index_file):
    cdef:
        str line, chrom, strand, key
        int start, end
        float factor

    with open(index_file,'r') as idx:
        index = pickle.load(idx)


    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end ,strand = parse_bed_line(line)
            seq = genome_fa.fetch(chrom, start, end)
            seq = reverse_complement(seq) if strand == "-" else seq
            key = seq[:3] + ','+ seq[-3:]
            factor = index[key.upper()]
            print(line.strip() + '\t%.3f' %(factor), file=sys.stdout)