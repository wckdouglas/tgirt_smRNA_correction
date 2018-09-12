
import pysam
import numpy as np
import pandas as pd
from operator import itemgetter
from sequencing_tools.fastq_tools import reverse_complement
from collections import defaultdict
from functools import reduce
from itertools import product
import pickle

def merge_trinucleotide(x, y):
    return x.merge(y, on ='trinucleotide', how='outer') 

   

class build_weights():
    def __init__(self, args):
        '''
        take a bed file, and a fasta file:

        # follow hexamer correction but only use 3 nucleotides
        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2896536/
        bed_file = '/stor/work/Lambowitz/cdw2854/miRNA/ercc/ERCC.bed'
        fa_file = '/stor/work/Lambowitz/ref/RNASeqConsortium/ercc/ercc.fa'
        '''

        self.bed_file = args.bed
        self.fa = pysam.Fastafile(args.fasta)
        self.base_dict = None
        self.base_df = None


    def analyze_bed_ends(self, max_iter = 100000):
        cdef:
            int i = 0
            str line, chrom, start, end, strand
            str seq, subseq
        
        self.base_dict = defaultdict(lambda: defaultdict(int))
        with open(self.bed_file) as inbed:
            while i < max_iter:
                try:
                    line = next(inbed) 
                    fields = line.strip().split('\t')
                    chrom,start,end,strand = itemgetter(0,1,2,5)(fields)
                    seq = self.fa.fetch(chrom, int(start), int(end))
                    seq = reverse_complement(seq) if strand == '-' else seq
                    
                    
                    self.base_dict['5'][seq[:3]] += 1
                    self.base_dict['3'][seq[-3:]] += 1

                    subseq = seq[3:-3]
                    for (b1,b2,b3) in zip(subseq, subseq[1:], subseq[2:]):
                        self.base_dict['background'][b1+b2+b3] += 1
                    i += 1
                except StopIteration:
                    break
            

        self.base_df = pd.DataFrame()\
            .from_dict(self.base_dict)\
            .transform(lambda x: x/x.sum(axis=0))\
            .reset_index() \
            .assign(w3 = lambda d: d['3']/d['background'])\
            .assign(w5 = lambda d: d['5']/d['background']) 


    def base_dict_to_weights(self):
        self.weight_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))

        for i, row in self.base_df.iterrows():
            self.weight_dict['3'][row['index']] = row['w3']
            self.weight_dict['5'][row['index']] = row['w5']


    def output_weights(self):
        cdef:
            str tail, head

        self.index = {}
        combination = [''.join(x) for x in product('ACTG',repeat=3)]
        for tail in combination:
            tail_score = self.weight_dict['3'][tail]
            for head in combination:
                head_score = self.weight_dict['5'][head]
                self.index[head + ',' + tail] = 1/np.sqrt(tail_score * head_score)  #geometric mean
    