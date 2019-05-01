from __future__ import print_function
import pickle
import sys
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
from itertools import product
from libc.math cimport exp, sqrt
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


def extract_kmer(str sequence, int k):
    '''
    output kmer
    '''
    cdef int i

    for i in range(len(sequence) - k + 1):
        yield(sequence[i:i+k])


class BuildWeights():
    def __init__(self, args):
        '''
        take a bam file:

        # follow hexamer correction but only use 3 nucleotides
        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896536/

        Usage:
        bw = build_weights(args)
        bw.analyze_bam_ends()
        bw.compute_weights()
        bw.output_weights()

        '''

        self.bam_file = args.inbam
        self.nucl = args.nucl
        self.background_start = self.nucl + 1
        self.weights_index = args.weight_index
        self.max_iter = args.iter
        self.base_df = None
        self.index = {}
        self.base_dict = defaultdict(lambda: defaultdict(int))
        self.weight_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))


    def analyze_bam_ends(self):
        '''
        analyze bam_reads
        First 3 nucleotides from read1 and read2 are store in different dictionary
        background is collected by looking at trinucloetides beyond position 3
        '''
        cdef:
            int read1_count = 0
            int read2_count = 0
            str end, strand, seq, subseq
            AlignedSegment aln
            str kmer
        
        pbar = tqdm(total=self.max_iter)
        print('Analyzing %i read pairs' %self.max_iter, file = sys.stderr)
        with pysam.Samfile(self.bam_file) as inbam:
            while read1_count < self.max_iter and read2_count < self.max_iter:
                try:
                    aln = next(inbam) 
                    seq = aln.get_forward_sequence()
                    
                    if aln.is_read1:
                        read1_count += 1
                        end = 'read1'
                    else: 
                        end = 'read2'
                        read2_count += 1
                        pbar.update(1)

                    self.base_dict[end][seq[:self.nucl]] += 1

                    subseq = seq[self.background_start:]
                    for kmer in extract_kmer(subseq, self.nucl):
                        self.base_dict['background_' + end][kmer] += 1
                except StopIteration:
                    break

        self.__base_dict_to_weights__()
            

    def __base_dict_to_weights__(self):
        ## scores are in log scale odd
        self.base_df = pd.DataFrame()\
            .from_dict(self.base_dict)\
            .transform(lambda x: np.log(x) - np.log(x.sum(axis=0)))\
            .reset_index() \
            .assign(read1_weights = lambda d: d['background_read1'] - d['read1'],
                    read2_weights = lambda d: d['background_read2'] - d['read2']) 

        for i, row in self.base_df.iterrows():
            self.weight_dict['read1'][row['index']] = row['read1_weights']
            self.weight_dict['read2'][row['index']] = row['read2_weights']


    def compute_weights(self):
        '''
        indexing the weights,
        it is the geometric mean between read1 weight and read2 weight
        everything is in log odd before putting into index
        '''
        cdef:
            str read1_seq, read2_seq
            double read1_weight, read2_weight

        combination = [''.join(x) for x in product('ACTG',repeat=self.nucl)]
        for read1_seq in combination:
            read1_weight = self.weight_dict['read1'][read1_seq]
            for read2_seq in combination:
                read2_weight = self.weight_dict['read2'][read2_seq]
                self.index[read1_seq + ',' + read2_seq] = exp(0.5 * (read1_weight + read2_weight))  

    def output_weights(self):
        '''
        save index
        '''
        with open(self.weights_index, 'wb') as index:
            pickle.dump(self.index, index)
        print('Written %s' %self.weights_index, file = sys.stderr)
