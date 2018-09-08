#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import explained_variance_score
import os
import sys
from itertools import product
import pickle
from scipy.stats import pearsonr

def positioning(x):
    return x[-1]

def count_to_cpm(count_array):
    count_array = np.true_divide(count_array,count_array.sum()) * 1e6 
    return count_array

def get_end(x):
    if 'head' in x:
        return "5'"
    elif 'tail' in x:
        return "3'"

def make_column_name(colnames):
    col_d = pd.DataFrame({'nucleotide':colnames.str.slice(-1),
             'position':colnames.str.slice(4,5),
             'end':colnames.map(get_end)}) \
        .assign(offset = lambda d: np.where(d.end=="5'",-1, 3)) \
        .assign(adjusted_position = lambda d: np.abs(d.position.astype(int) - d.offset))\
        .assign(colnames = colnames)
    return col_d.end + '-position:'+col_d.adjusted_position.astype(str) +':'+ col_d.nucleotide

def preprocess_dataframe(df):
    nucleotides = df.columns[df.columns.str.contains('head|tail')]
    dummies = pd.get_dummies(df[nucleotides])
    dummies.columns = make_column_name(dummies.columns)
    df = pd.concat([df,dummies],axis=1) \
        .drop(nucleotides, axis=1) 
    return df

class lm_model():
    '''

    Train a ridge model for predicting log2 cpm deviation from golden standard
    '''
    def __init__(self, train_set, index_file, num_nucleotide):
        self.index = {}
        self.train_set = train_set
        self.index_file = index_file
        self.lm = Ridge()
        self.num_nucleotide


    def preprocess_data(self):
        '''
        make dummy variable for end nucleotide, and
        indicate training set:

        X: nx24 matrix, every positional nucleotide as columns
        Y: log2 CPM difference: experimetnal count - expected count  
        '''
        self.df = pd.read_table(self.train_set)\
            .assign(expected_cpm = lambda d: count_to_cpm(d['expected_count']))\
            .assign(experimental_cpm = lambda d: count_to_cpm(d['experimental_count']))\
            .assign(log_cpm_diff = lambda d: np.log2(d.experimental_cpm+1) - np.log2(d.expected_cpm+1))\
            .pipe(preprocess_dataframe) \
            .drop(['experimental_count','experimental_cpm','expected_cpm','expected_count'], axis=1)

        self.X = self.df.drop(['seq_id','log_cpm_diff'], axis=1)
        self.Y = self.df['log_cpm_diff'].values
    

    def train_lm(self):
        '''
        Train a ridge model for correction
        '''
        self.lm.fit(self.X, self.Y)
        pred_Y = self.lm.predict(self.X)
        rsqrd = explained_variance_score(self.Y, pred_Y)
        rho, pval = pearsonr(pred_Y, self.Y)
        print('Trained model', file=sys.stderr)
        print('R-sqrd: %.2f' %(rsqrd), file=sys.stderr)
        print('Pearson correlation: %.2f' %(rho), file=sys.stderr)


    def write_index(self):
        '''

        For each combination of 3' and 5' trinucleotide, generate a correction factor

        '''
        coef_dict = {n:c for c, n in zip(self.lm.coef_,self.X.columns)}
        combination = [''.join(x) for x in product('ACTG',repeat=self.num_nucleotide)]

        for tail in combination:
            tail_score = sum(coef_dict["3'-position:%i:%s" %(self.num_nucleotide-i,t)] for i, t in enumerate(tail))
            for head in combination:
                head_score = sum(coef_dict["5'-position:%i:%s" %(i+1,t)] for i, t in enumerate(head))
                self.index[head + ',' + tail] = tail_score + head_score
        
        with open(self.index_file,'wb') as idx_file:
            pickle.dump(self.index, idx_file)


        print('Make index: %s' %self.index_file, file=sys.stdout)
