from __future__ import print_function
import pickle
import sys
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
from itertools import product
from libc.math cimport exp
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


def merge_trinucleotide(x, y):
    return x.merge(y, on ='trinucleotide', how='outer') 


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
        self.weights_index = args.weight_index
        self.max_iter = args.iter
        self.base_df = None
        self.index = {}
        self.base_dict = defaultdict(lambda: defaultdict(int))
        self.weight_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))


    def analyze_bam_ends(self):
        '''
        analyze bam_reads
        '''
        cdef:
            int read1_count = 0
            int read2_count = 0
            str line, chrom, start, end, strand
            str seq, subseq
            AlignedSegment aln
        
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

                    self.base_dict[end][seq[:3]] += 1

                    subseq = seq[5:]
                    for (b1,b2,b3) in zip(subseq, subseq[1:], subseq[2:]):
                        self.base_dict['background'][b1+b2+b3] += 1
                except StopIteration:
                    break
            

        ## scores are in log scale odd
        self.base_df = pd.DataFrame()\
            .from_dict(self.base_dict)\
            .transform(lambda x: np.log(x) - np.log(x.sum(axis=0)))\
            .reset_index() \
            .assign(read1_weights = lambda d: d['background'] - d['read1'])\
            .assign(read2_weights = lambda d: d['background'] - d['read2']) 


    def base_dict_to_weights(self):

        for i, row in self.base_df.iterrows():
            self.weight_dict['read1'][row['index']] = row['read1_weights']
            self.weight_dict['read2'][row['index']] = row['read2_weights']


    def compute_weights(self):
        cdef:
            str tail, head

        combination = [''.join(x) for x in product('ACTG',repeat=3)]
        for tail in combination:
            tail_score = self.weight_dict['read1'][tail]
            for head in combination:
                head_score = self.weight_dict['read2'][head]
                self.index[head + ',' + tail] = exp(0.5 * (tail_score + head_score) )  #geometric mean

    def output_weights(self):
        with open(self.weights_index, 'wb') as index:
            pickle.dump(self.index, index)
        print('Written %s' %self.weights_index, file = sys.stderr)