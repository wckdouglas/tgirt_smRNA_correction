#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import GridSearchCV
import os
import sys
from itertools import product
import pickle
from scipy.stats import pearsonr
from collections import defaultdict
import re
from .model import h2o_rf

def count_to_cpm(count_array):
    count_array = np.true_divide(count_array,count_array.sum()) * 1e6 
    return count_array



class nucleotide_model():
    '''

    Train a ridge model for predicting log2 cpm deviation from golden standard
    '''
    def __init__(self, train_set, index_file, num_nucleotide):
        self.index = {}
        self.train_set = train_set
        self.index_file = index_file
        self.coef_file = index_file.replace('_index.pkl','_coef.pkl')
        self.model = h2o_rf
        self.num_nucleotide = num_nucleotide


    def preprocess_data(self):
        '''
        make dummy variable for end nucleotide, and
        indicate training set:

        X: nx24 matrix, every positional nucleotide as columns
        Y: log2 CPM difference: experimetnal count - expected count  
        '''
        self.train_df = pd.read_table(self.train_set)\
            .assign(expected_cpm = lambda d: count_to_cpm(d['expected_count']))\
            .assign(experimental_cpm = lambda d: count_to_cpm(d['experimental_count']))\
            .query('experimental_cpm > 0')\
            .assign(log_cpm_diff = lambda d: np.log(d.experimental_cpm) - np.log(d.expected_cpm))\
            .drop(['experimental_count','experimental_cpm','expected_cpm','expected_count'], axis=1)

        #self.X = self.df.drop(['seq_id','log_cpm_diff'], axis=1)
        self.X = self.df.filter(regex="^head|^tail")
        self.Y = self.df['log_cpm_diff'].values


    def train_rf(self):
        '''
        Train a ridge model for correction
        '''
        self.model.fit(self.X, self.Y)
        pred_Y = self.model.predict(self.X)
        rsqrd = r2_score(self.Y, pred_Y)
        rho, pval = pearsonr(pred_Y, self.Y)
        print('Trained model', file=sys.stderr)
        print('R-sqrd: %.2f' %(rsqrd), file=sys.stderr)
        print('Pearson correlation: %.2f' %(rho), file=sys.stderr)

   

    def train_lm(self):
        '''
        Train a ridge model for correction
        '''
        self.lm.fit(self.X, self.Y)
        pred_Y = self.lm.predict(self.X)
        rsqrd = r2_score(self.Y, pred_Y)
        rho, pval = pearsonr(pred_Y, self.Y)
        print('Trained model', file=sys.stderr)
        print('R-sqrd: %.2f' %(rsqrd), file=sys.stderr)
        print('Pearson correlation: %.2f' %(rho), file=sys.stderr)


    def write_index(self):
        '''

        For each combination of 3' and 5' trinucleotide, generate a correction factor

        '''
        with open(self.index_file, 'wb') as index_file:
            model_param = {'model': self.model, 
                        'X_col': self.X.columns.tolist()}
            pickle.dump(model_param, index_file)

        print('Make index: %s' %self.index_file, file=sys.stdout)
