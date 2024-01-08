#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import itertools
import HLAGuessr.model_tools as mt
import os


class PreprocessedModel():
    
    def __init__(self,chain,sig_threshold=0.05, seq_threshold=None, list_file=None, hla_file=None, files_path=None,weights_file=None):
        
        self.main_directory = os.path.dirname(__file__)
        
        if any([x is not None for x in [list_file, hla_file, files_path,weights_file]]):
            self.list_file, self.hla_file, self.files_path, self.weights_file = list_file, hla_file, files_path, weights_file
        else:
            self.list_file = self.main_directory+'/Training_data/inference_all_data_TCRs_sign_HLA.tsv'
            self.hla_file = self.main_directory+'/Training_data/HLA_grouped_patients.tsv'
            self.weights_file = self.main_directory+'/Training_data/HLA_frequency_for_validation.tsv'

        self.chain = chain
        
        if len(self.chain)==2:
                self.reg_file_path = self.main_directory + '/Training_data/classifier_params_alpha+beta.tsv'
        if len(self.chain)==1:
            if 'alpha' in self.chain:
                self.reg_file_path = self.main_directory + '/Training_data/classifier_params_alpha+beta.tsv'
            if 'beta' in self.chain:
                self.reg_file_path = self.main_directory + '/Training_data/classifier_params_beta.tsv'
        self.df_param = pd.read_csv(self.reg_file_path,delimiter='\t')
            
        self.data_train = self.load_training_data()
        self.t, self.n_seq = sig_threshold, seq_threshold
        
    def load_training_data(self):
        
        big_df = pd.DataFrame()
        for c in self.chain:
            f = self.main_directory+'/Training_data/{}_cdr3aa_shared_at_least_3.txt'.format(c)
            df = pd.read_csv(f, delimiter='\t')
            df['chain'] = c
            big_df = pd.concat([big_df, df], ignore_index=True)
        big_df.set_index('cdr3+v_family',inplace=True)
        return(big_df)
    
    def get_hla_info(self,hla_threshold,hla_target):
        df = pd.read_csv(self.hla_file,delimiter='\t')
        df_hla = df.loc[((df['n_positive']>hla_threshold)&(df['HLA']==hla_target)),:]
        return(df_hla)
    
    def train_patients_database(self, hla_target, df_hla):  #,untyped,pred_files
        
        if 'beta' not in self.chain:
            p_total = [p for p in df_hla[df_hla['HLA']==hla_target].pos_patients.values[0].split(', ')+df_hla[df_hla['HLA']==hla_target].neg_patients.values[0].split(', ') if '_Em' not in p]
        else:
            p_total = df_hla[df_hla['HLA']==hla_target].pos_patients.values[0].split(', ')+df_hla[df_hla['HLA']==hla_target].neg_patients.values[0].split(', ')

        n_patients = len(p_total)
        values = list(np.arange(n_patients))
        dic = {p_total[i]: values[i] for i in range(len(p_total))}
            
        return(dic,n_patients)

    def get_significant_tcr(self, hla_target):
        
        if len(self.chain)==1:
            df = pd.read_csv(self.list_file, delimiter='\t')
            df = df.loc[((df['p_BH']<self.t)&(df['chain']==self.chain[0])&(df['Attribute']=='Overrepresented')), ['TCR+V','p_BH','HLA']]
        if len(self.chain)==2:
            df = pd.read_csv(self.list_file, delimiter='\t')
            df_alpha = df.loc[((df['p_BH']<self.t)&(df['chain']==self.chain[0])&(df['Attribute']=='Overrepresented')), ['TCR+V','p_BH','HLA']]
            df = pd.read_csv(self.list_file, delimiter='\t')
            df_beta = df.loc[((df['p_BH']<self.t)&(df['chain']==self.chain[1])&(df['Attribute']=='Overrepresented')), ['TCR+V','p_BH','HLA']]
            if self.n_seq:
                df_alpha = df_alpha[:self.n_seq_train]
                df_beta = df_beta[:self.n_seq_train]
            df = pd.concat([df_alpha,df_beta],ignore_index=True)
        sign_tcr = df[df['HLA']==hla_target]['TCR+V'].tolist()
        return(sign_tcr)

    def look_for_tcr_index(self,df,seqs):
        
        idx1 = df.index
        idx2 = pd.Index(seqs)
        intersection = idx1.intersection(idx2)
        df1 = df.loc[intersection]
        df1['cdr3+v_family'] = df1.index
        df1.reset_index(drop=True,inplace=True)
        return(df1)

    def get_occurrence_matrix_train(self,data,n_patients,dic,grouped,ref_OMs=None,matrix_index=None):
           
        OMs = []
        data.drop_duplicates('cdr3+v_family',keep='first',inplace=True)
        data.sort_values('cdr3+v_family',inplace=True)
        for i,c in enumerate(self.chain):
            df1 = data[(data['chain']==c)].reset_index()
            if grouped==True:
                m = int(n_patients) 
                n = int(len(df1))
                OM = np.zeros([m,n])
                for i,x in enumerate(df1['Patient']):
                    for y in x.split(', '):
                        if y in dic:
                            j = dic.get(str(y))
                            OM[int(j)][i] = 1
                OMs.append(OM)
        
        z = None
        if len(self.chain)==2:
            z = np.zeros([m,1])
            for j,k in enumerate(dic): # Add z parameter as a feature
                if '_Em' in k:
                    z[j][0]=1
        return(OMs,z)
    
    def get_occurrence_matrix_test(self,OMs_train,n_patients,matrix_index):
        
        OMs = []
        for i,c in enumerate(self.chain):
            m = int(n_patients) 
            n = int(len(OMs_train[i][0]))
            OM = np.zeros([m,n])
            matrix_index_c = matrix_index[matrix_index['chain']==c]
            for k in list(matrix_index_c.index):
                OM[0][k]=1
            OMs.append(OM)
        z = np.zeros([m,1])
        return(OMs,z)
    
    def get_matches_x(self,M,z=None):
        
        if 'alpha' in self.chain and 'beta' in self.chain:
            X = np.concatenate((M[0],M[1],z), axis=1)
        if 'alpha' not in self.chain or 'beta' not in self.chain:
            X = M[0]
        return(X)
    
    def get_matches_y(self, n_patients, dic, df_hla):
        
        OH = np.zeros([1,int(n_patients)])
        for i,x in enumerate(df_hla['pos_patients'].values):
            for c in str(x).split(', '):
                if c in dic:
                    if 'beta' in self.chain:
                        i = dic.get(str(c)) 
                    else:
                        if '_Em' not in c:
                            i = dic.get(str(c))
                    OH[0][int(i)]=1
        y = list(OH[0])
        return(y)