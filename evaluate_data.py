#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import HLAGuessr.model_tools as mt
from sklearn.model_selection import train_test_split,ParameterGrid
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import roc_curve,roc_auc_score,RocCurveDisplay,ConfusionMatrixDisplay,classification_report,confusion_matrix,accuracy_score
from HLAGuessr.load_model import PreprocessedModel

class GetProbabilities(object):  # I need to modify the input format to cdr3aa+v_family
    
    def __init__(self,hla_target,chain,hla_threshold=2,data_test=None):
        
        self.chain = chain
        self.data_test = data_test
        pm = PreprocessedModel(self.chain)
        self.df_hla = pm.get_hla_info(hla_threshold,hla_target)
        self.dic_p_train,self.n_pat_train = pm.train_patients_database(hla_target,self.df_hla)
        self.dic_p_test,self.n_pat_test = self.test_patients_database(hla_target, self.data_test)
        self.sign_tcr = pm.get_significant_tcr(hla_target)
        self.df_w_sign_train = pm.look_for_tcr_index(pm.data_train, self.sign_tcr)
        self.df_w_sign_test = pm.look_for_tcr_index(self.data_test, self.sign_tcr)
        if len(self.df_w_sign_test)==0:
            print('No HLA-related TCRs have been found in your repertoire')
            print('Exiting...')
            break
        else:
            if len(self.chain)==2:
                self.OMs_train, self.z_train = pm.get_occurrence_matrix(self.df_w_sign_train,self.n_pat_train,self.dic_p_train,grouped=True)
                self.OMs_test, self.z_test = pm.get_occurrence_matrix(self.df_w_sign_test,self.n_pat_test,self.dic_p_test,grouped=False)
                self.X_train = pm.get_matches_x(self.OMs_train,self.z_train)
                self.X_test = pm.get_matches_x(self.OMs_test,self.z_test)
                if len(self.chain)==1:
                    self.OMs_train = pm.get_occurrence_matrix(self.df_w_sign_train,self.n_pat_train,self.dic_p_train,grouped=True)
                    self.OMs_test = pm.get_occurrence_matrix(self.df_w_sign_test,self.n_pat_test,self.dic_p_test,grouped=False)
                    self.X_train = pm.get_matches_x(self.OMs_train)
                    self.X_test = pm.get_matches_x(self.OMs_test) 
                self.y_train = pm.get_matches_y(self.n_pat_train,self.dic_p_train,self.df_hla)
                self.X_scaled = self.renormalize_training_matrix(hla_target,self.df_w_sign_train,self.df_w_sign_test,self.X_train)
                self.hla_prob = self.get_hla_probabilities(self.X_scaled,self.y_train,self.X_test,hla_target,pm.df_param)
    
    def test_patients_database(self, hla_target, data_test):  #,untyped,pred_files
        p = list(set(list(data_test['Patient'])))
        n_patients = len(p)
        dic = {p[i]: list(np.arange(n_patients))[i] for i in range(len(p))}
        return(dic,n_patients)
        
    def renormalize_training_matrix(self,hla_target, df_w_sign_train, df_w_sign_test,X_train):
        
        df_w_sign_train.reset_index(drop=True,inplace=True)
        df_w_sign_test.drop_duplicates('cdr3+v_family',keep='first',inplace=True)
        df_w_sign_test.sort_values('cdr3+v_family',inplace=True)
        sign_tcr_2 = list(df_w_sign_test['cdr3+v_family'])
        matrix_index = mt.applyParallel(sign_tcr_2, df_w_sign_train, mt.get_tcr_index)
        l = list(matrix_index.index)
        l.extend([len(X_train[0])-1])
        X_norm = X_train[:,l]
        return(X_norm)
    
    def get_hla_probabilities(self,X_train_scaled,y_train,X_test,hla_target,df_param,solver='lbfgs',penalty='l2',l1_ratio=None,class_weight=None,c=None):
        
        if c is None:
            c = np.float(df_param[df_param['HLA']==hla_target]['Regularization'])
        
        model_lr = LogisticRegression(solver=solver,penalty=penalty,C=c)
        model_lr.fit(X_train_scaled, y_train)
        prob_lr = model_lr.predict_proba(X_test)
        return(prob_lr[0][1])
    
    def classifier_params(self,df_param,hla_target):
        df1 = df_param[df_param['HLA']==hla_target]
        auc = np.float(df1['AUC'])
        a = np.float(df1['accuracy'])
        p = np.float(df1['precision'])
        s1 = np.float(df1['sensitivity'])
        s2 = np.float(df1['specificity'])
        return([auc,a,p,s1,s2])