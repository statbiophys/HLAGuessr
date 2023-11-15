#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


translate_header = {
    'aminoAcid': 'cdr3aa',
    'cdr3_amino_acid': 'amino_acid',
    'vGeneName': 'v_gene',
    'junction_aa': 'cdr3aa',
    'CDR3.amino.acid.sequence':'cdr3aa',
    'amino_acid':'cdr3aa',
    'cdr3_aa':'cdr3aa',
    'aa':'cdr3aa',
    'bestVGene':'v_gene',
    'v':'v_gene',
    'v_family':'v_gene'}

class Processing(object):
    
    def __init__(self,chain,alpha_input=None,beta_input=None,delimiter='\t'):

        #self.read_thresh=read_thresh
        #self.max_length=max_length
        self.alpha_input = alpha_input
        self.beta_input = beta_input
        self.delimiter = delimiter
        self.chain=chain
        self.delimiter = delimiter
        self.data_test = self.load_test_data(self.alpha_input,self.beta_input,self.delimiter)
        
    def load_test_data(self,alpha_files,beta_files,sep):

        alpha = pd.DataFrame()
        beta = pd.DataFrame()

        if alpha_files is not None:
            alpha = pd.read_csv(alpha_files,delimiter=sep)
            format_alpha = self.format_dataframe(alpha)
            format_alpha['chain'] = 'alpha'
            
        if beta_files is not None:
            beta = pd.read_csv(beta_files,delimiter='\t')
            format_beta = self.format_dataframe(beta)
            format_beta['chain'] = 'beta'

        big_df = pd.concat([format_alpha,format_beta], ignore_index=True)
        big_df.set_index('cdr3+v_family',inplace=True)
        return(big_df)
    
    def format_dataframe(self,df):
        '''
        Implement filtering pipeline
        '''
        new_columns = self.rename_cols(df)
        df.columns = new_columns
        new_v_gene = self.format_genes(df)
        df['v_gene'] = new_v_gene
        
        prod_df = self.select_productive(df)
        #self.select_cdr3_length()
        prod_df['cdr3+v_family'] = prod_df['cdr3aa']+'+'+prod_df['v_gene']
        prod_df = prod_df.loc[:,['cdr3+v_family','Patient']]

        return(prod_df)
            
    def format_genes(self,df):   
        '''
        Considers only v gene family
        '''
        new_vs=df.v_gene.apply(lambda x: x.split('*')[0]).values
        return(list(new_vs))
        
    def select_productive(self,df):
        '''
        looks for *,_,nan,~ in the cdr3.
        '''
        bad_cdr3=df.cdr3aa.str.contains('\*|_|~',na=True,regex=True)
        df = df[~bad_cdr3]
        df.reset_index(inplace=True,drop=True)
        return(df)
        
    def rename_cols(self,df):
        columns=df.columns.values
        new_cols=[]
        for col in columns:
            try: new_cols.append(translate_header[col])
            except: new_cols.append(col)
        return(new_cols)
    
    def select_cdr3_length(self):
        '''
        select cdr3 length smaller than max_length
        '''
        self.df['selection_length']=self.df.amino_acid.fillna('nan').apply(len)<self.max_length
        print ('long cdr3s:',np.sum(np.logical_not(self.df['selection_length'].values[self.df['selection'].values])))
        self.df['selection']=np.logical_and(self.df['selection'].values,self.df['selection_length'].values)
