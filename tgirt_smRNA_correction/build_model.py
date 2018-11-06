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

def make_column_name(colnames, num_nucleotide):
    col_d = pd.DataFrame({'nucleotide':colnames.str.slice(-1),
             'position':colnames.str.slice(4,5),
             'end':colnames.map(get_end)}) \
        .assign(offset = lambda d: np.where(d.end=="5'",-1, num_nucleotide)) \
        .assign(sign = lambda d: np.where(d.end=="5'",'+','-')) \
        .assign(adjusted_position = lambda d: np.abs(d.position.astype(int) - d.offset))\
        .assign(colnames = colnames)
    return col_d.end + '-position:N' + col_d.sign +col_d.adjusted_position.astype(str) +':'+ col_d.nucleotide



def preprocess_dataframe(df, num_nucleotide):
    nucleotides = df.columns[df.columns.str.contains('^head[0-9]+$|^tail[0-9]+')]
    dummies = pd.get_dummies(df[nucleotides])
    dummies.columns = make_column_name(dummies.columns, num_nucleotide)
    cols = dummies.columns.tolist()
    cols.sort()
    dummies = dummies.loc[:, cols]

    df = pd.concat([df,dummies],axis=1) \
        .drop(nucleotides, axis=1) 
    return df


def preprocess_rf_dataframe(df, num_nucleotide):
    '''
    labelEncoder for categorical training in RF
    '''
    le = LabelEncoder()
    le.fit(list('ACTG'))
    nucleotides = df.columns[df.columns.str.contains('^head[0-9]+$|^tail[0-9]+$')]
    for col in nucleotides:
        pos_end = "5'" if 'head' in col else "3'"
        pos = int(col[-1])
        pos = pos + 1 if pos_end == "5'" else pos - num_nucleotide
        new_col = pos_end +':'+ str(pos)
        df[new_col] = le.transform(df[col])
    
    return df


class nucleotide_model():
    '''

    Train a ridge model for predicting log2 cpm deviation from golden standard
    '''
    def __init__(self, train_set, index_file, num_nucleotide):
        self.index = {}
        self.train_set = train_set
        self.index_file = index_file
        self.coef_file = index_file.replace('_index.pkl','_coef.pkl')
        rf = RandomForestRegressor()

        mdepth = list(range(90,100))
        mdepth.append(None)
        self.rf = GridSearchCV(estimator=rf,  
                cv = 8, n_jobs=-1,
                param_grid={'min_samples_split':range(2,10),
                            'max_depth' : mdepth,
                            'n_estimators': range(10,20)}) 
        self.num_nucleotide = num_nucleotide


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
            .query('experimental_cpm > 0')\
            .assign(log_cpm_diff = lambda d: np.log(d.experimental_cpm) - np.log(d.expected_cpm))\
            .pipe(lambda d: preprocess_rf_dataframe(d, self.num_nucleotide)) \
            .drop(['experimental_count','experimental_cpm','expected_cpm','expected_count'], axis=1)

        #self.X = self.df.drop(['seq_id','log_cpm_diff'], axis=1)
        self.X = self.df.filter(regex="^[35]':")
        self.Y = self.df['log_cpm_diff'].values


    def train_rf(self):
        '''
        Train a ridge model for correction
        '''
        self.rf.fit(self.X, self.Y)
        pred_Y = self.rf.predict(self.X)
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
            model_param = {'model': self.rf, 
                        'X_col': self.X.columns.tolist()}
            pickle.dump(model_param, index_file)

        print('Make index: %s' %self.index_file, file=sys.stdout)
