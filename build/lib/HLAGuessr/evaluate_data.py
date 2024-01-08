#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import HLAGuessr.model_tools as mt
from sklearn.model_selection import train_test_split,ParameterGrid, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import (roc_curve,roc_auc_score,RocCurveDisplay,ConfusionMatrixDisplay,classification_report,confusion_matrix,average_precision_score,log_loss,f1_score)
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
            self.hla_prob = 'N/A'
            self.hla_output = False
        else:
            self.OMs_train, self.z_train = pm.get_occurrence_matrix_train(self.df_w_sign_train,self.n_pat_train,self.dic_p_train,True)
            self.X_train = pm.get_matches_x(self.OMs_train,self.z_train)
            self.X_scaled, self.matrix_index = self.renormalize_training_matrix(hla_target,self.df_w_sign_train,self.df_w_sign_test,self.X_train,self.chain)
            self.OMs_test, self.z_test = pm.get_occurrence_matrix_test(self.OMs_train,self.n_pat_test,self.matrix_index)
            self.X_test = pm.get_matches_x(self.OMs_test, self.z_test)
            self.y_train = pm.get_matches_y(self.n_pat_train,self.dic_p_train,self.df_hla)
            self.class_weights = mt.get_weights(pm.weights_file,hla_target)
            self.hla_prob, self.hla_output = self.get_hla_probabilities(self.X_train,np.array(self.y_train),self.X_test,hla_target,self.class_weights)
    
    def test_patients_database(self, hla_target, data_test):  #,untyped,pred_files
        p = list(set(list(data_test['Patient'])))
        n_patients = len(p)
        dic = {p[i]: list(np.arange(n_patients))[i] for i in range(len(p))}
        return(dic,n_patients)
        
    def renormalize_training_matrix(self,hla_target, df_w_sign_train, df_w_sign_test,X_train,chain):
        l = []
        big_matrix_index = pd.DataFrame()
        for c in chain:
            df_w_sign_train_c = df_w_sign_train.loc[df_w_sign_train['chain']==c,:]
            df_w_sign_train_c.reset_index(drop=True,inplace=True)
            df_w_sign_test.drop_duplicates('cdr3+v_family',keep='first',inplace=True)
            df_w_sign_test.sort_values('cdr3+v_family',inplace=True)
            sign_tcr_2 = list(df_w_sign_test['cdr3+v_family'])
            matrix_index = mt.applyParallel(sign_tcr_2, df_w_sign_train_c, mt.get_tcr_index)
            l0 = list(matrix_index.index)
            big_matrix_index = pd.concat([big_matrix_index,matrix_index])
            l = l+l0
            
        if 'alpha' in chain and 'beta' in chain:
            l.extend([len(X_train[0])-1])
        X_norm = X_train[:,l]
        return(X_norm,big_matrix_index)
    
    def get_hla_probabilities(self,X_train_scaled, y_train,X_test,hla_target,class_weight,solver='liblinear',penalty='l1',t=0.9005,l1_ratio=None):
        paramGrid = {
    'C': [50,20,10]+ list(np.logspace(-5,0,50)),
    'solver':[solver], 
    'penalty': [penalty],
    'class_weight':[class_weight]
        }
        lr = LogisticRegression()
        model_lr = GridSearchCV(lr, param_grid = paramGrid, cv = 5, verbose=True, n_jobs=-1, scoring='f1_macro') #neg_log_loss
        
        model_lr.fit(X_train_scaled, y_train)
        prob_lr = model_lr.predict_proba(X_test)
        predictions_lr = (prob_lr[0][1] > t) #.astype(int)
        return(prob_lr[0][1],predictions_lr)
    
    def classifier_params(self,df_param,hla_target):
        df1 = df_param[df_param['HLA']==hla_target]
        auc = float(df1['AUC'])
        a = float(df1['accuracy'])
        p = float(df1['precision'])
        s1 = float(df1['sensitivity'])
        s2 = float(df1['specificity'])
        return([auc,a,p,s1,s2])