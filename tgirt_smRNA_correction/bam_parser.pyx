from __future__ import print_function
import pysam
import pickle
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


def read_paired_bam(AlignmentFile bam_handle):
    '''
    paired end bam file generator
    '''
    cdef:
        AlignedSegment read1
        AlignedSegment read2

    try:
        while True:
            #try:
                read1 = next(bam_handle)
                read2 = next(bam_handle)

                assert read1.query_name.split('/')[0] == read2.query_name.split('/')[0], 'Read1 and Read2 not next to each other?'
                read1,read2 = (read2,read1) if read2.is_read1 else (read1,read2)
                assert(read1.is_read1 and read2.is_read2)
                yield read1, read2

    except StopIteration:
        pass


class BAMCorrector():
    '''
    adding an "ZW" tag for weights for each bam alignment
    '''
    def __init__(self, args):
        self.bam_file = args.inf
        self.index = args.index
        self.out_bam = args.outf
        self.tag = args.tag

        with open(self.index, 'rb') as weights:
            self.index_table = pickle.load(weights)


    def correction(self):
        '''
        adding the tag to the alignments

        1. read both ends of a read pair
        2. collect the end sequence
        3. fetch weight from prebuilt weight index
        4. add to ZW tag
        '''

        cdef:
            AlignedSegment read1, read2
            str head , tail
            double weight
            int paired_count, thrown=0

        with pysam.Samfile(self.bam_file) as inbam:
            with pysam.Samfile(self.out_bam, 'w', template = inbam) as outbam:
                for paired_count, (read1, read2) in enumerate(read_paired_bam(inbam)):
                    head = read1.get_forward_sequence()[:3]
                    tail = read2.get_forward_sequence()[:3]
                    if 'N' not in head and 'N' not in tail:
                        weight = self.index_table[head + ',' + tail]
                        read1.set_tag('ZW', weight, value_type = 'f', replace=True)
                        read2.set_tag('ZW', weight, value_type = 'f', replace=True)
                        outbam.write(read1)
                        outbam.write(read2)
                    else:
                        thrown += 1

        print('Reweighted %i read pairs, threw away %i with "N" base' %(paired_count, thrown))


