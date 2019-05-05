#!/usr/bin/env python

from __future__ import print_function
import pysam
import argparse
import pyximport
import sys
import os
import pickle
import pkg_resources
model_dir = pkg_resources.resource_filename('tgirt_smRNA_correction', 'model')
default_model = model_dir + '/weights.pkl'


def get_opt():
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]))
    subparser = parser.add_subparsers(help='Reweighting TGIRT-seq reads using end sequence biases. subcommand:', dest = 'subcommand')
    subparser.required=True
        
    #add build index
    index_builder = subparser.add_parser('build', 
                    help='Building nucleotide bias table for each transcripts')
    index_builder.add_argument('-f','--fasta', help='Small RNA fasta file', required=True)
    index_builder.add_argument('-n','--nucleotide', help='How many nucleotide to look at from both end?', default=3, type=int)
    index_builder.add_argument('-o','--output_prefix', default='model',help='How many nucleotide to look at from both end?')
    index_builder.add_argument('-e','--expected_count', help='Expected count (comma delimintaed: seq,count) (default: all equal)')
    index_builder.add_argument('-c','--experimental_count', help='Experimental count (comma delimintaed: seq,count)', required=True)

    # weights builder
    train_weights =  subparser.add_parser('train', 
                    help='Building end nucleotide weights')
    train_weights.add_argument('-i', '--inbam', help='Input bam file', required=True)
    train_weights.add_argument('-x', '--weight_index', help = 'Output weight index (default: weight.pkl)', default='weight.pkl')
    train_weights.add_argument('-c', '--iter', help = 'How many reads to analyze for each end (default: 500000)', default=500000, type=int)
    train_weights.add_argument('-n', '--nucl', help = 'How many nucleotides from each end are used for computing the weight (default: 3)', default=3, type=int)
    train_weights.add_argument('--offset', help = 'How many nucleotides offset after biased-position (--nucl) to be considered as background (default: 0)', default=0, type=int)
    train_weights.add_argument('-l', '--len', help = 'How many nucleotides beyond --nucl is used for estimating background frequency (default: 50)', default=50, type=int)
    train_weights.add_argument('--debug', help = 'debug mode', action='store_true')
    
    # correction
    correction = subparser.add_parser('correct', 
            help='Using prebuilt bias index, add a "ZW" tag to bam alignments')
    correction.add_argument('-i', '--inf', help='Input fragment (name sorted bam file!!)', required=True)
    correction.add_argument('-x','--index', help = 'Weight index trained by "train" command (default: %s)' %default_model, default = default_model)
    correction.add_argument('-o','--outf', default='-',help='output file (default: -)')
    correction.add_argument('-t','--tag', default='ZW',help='BAM tag storing the weight (default: ZW)')
    correction.add_argument('--bed', action='store_true', 
                            help='input and output files are bed files'\
                                ', otherwise bam files')
    correction.add_argument('-f','--fasta', help='Genome fasta file, for correction with bed only')

    args = parser.parse_args()
    return args


def build(args):
    from tgirt_smRNA_correction.build_table import input_table
    from tgirt_smRNA_correction.build_model import nucleotide_model
    # building index    
    ## make nucleotide table
    base_table = input_table(args)
    base_table.build_base_table()
    base_table.add_observations()
    base_table.make_table()
    train_table = base_table.out_table_name


    ## train model
    index_file = args.output_prefix + '_index.pkl'
    model_file = args.output_prefix + '_model.pkl'
    model = nucleotide_model(train_table, index_file, args.nucleotide)
    model.preprocess_data()
    model.train_rf()
    model.write_index()


def train_weights(args):
    from tgirt_smRNA_correction.build_weights import BuildWeights
    bias_weights = BuildWeights(args)
    bias_weights.analyze_bam_ends()
    bias_weights.compute_weights()
    bias_weights.output_weights()


def correction(args):
    if args.bed:
        from tgirt_smRNA_correction.bed_parser import parse_bed, model_correction
        assert args.fasta is not None, 'Need a reference fasta file'
        fa = pysam.FastaFile(args.fasta)
        outfile = sys.stdout if args.outf == '-' or args.outf == '/dev/stdout' else open(args.outf,'w')

        idx = open(args.index,'rb')
        bias_index = pickle.load(idx) 
        idx.close()
        print('Retrieved index', file=sys.stderr)
        model_correction(args.inf, fa, bias_index, outfile)
        print('Added weights', file=sys.stderr)

    else:
        from tgirt_smRNA_correction.bam_parser import BAMCorrector
        assert args.inf.endswith('.bam'), 'Not a bam file?'
        bam_corrector = BAMCorrector(args)
        bam_corrector.correction()


def main():
    args = get_opt()

    if args.subcommand == 'build':
        build(args)
    
    elif args.subcommand == 'correct':
        correction(args)
    
    elif args.subcommand == 'train':
        train_weights(args)



if __name__ == '__main__':
    main()
