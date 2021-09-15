#!/usr/bin/env python

import os
import numpy as np
from collections import defaultdict
from tgirt_smRNA_correction.build_weights import extract_kmer, BuildWeights

class MakeTestCase():
    def __init__(self):
        self.test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'
        self.inbam = self.test_data_path + '/test.bam'
        self.nucl = 3
        self.offset = 4
        self.len = 20
        self.weight_index = self.test_data_path + '/weights.pkl'
        self.debug=False
        self.iter = 1000



def test_kmer():
    seq = 'ACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATA'
    kmer_counter = defaultdict(int)
    for i, kmer in enumerate(extract_kmer(seq, 3)):
        if i == 0:
            assert(kmer=="ACC")
        kmer_counter[kmer] += 1
    
    assert(kmer=="ATA")
    assert(kmer_counter['CAA'] == 5)


def test_weights():
    args = MakeTestCase()
    weight_builder = BuildWeights(args)
    weight_builder.analyze_bam_ends()
    assert(weight_builder.base_dict['read2']['GCC'] == 3)
    assert(weight_builder.base_dict['background_read2']['TGC'] == 63)
    assert(np.isclose(weight_builder.weight_dict['read1']['TGG'], -0.46851444))
