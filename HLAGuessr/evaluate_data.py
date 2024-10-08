#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd, numpy as np, HLAGuessr.model_tools as mt
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from HLAGuessr.load_model import PreprocessedModel
from HLAGuessr.processing import Processing

class GetProbabilities(object):
    
    def __init__(self,hla_target,chain,hla_threshold=None,ps=None,data_test=None,grouping=False,untyped=True,default_rep=True):
        if hla_threshold is None:
            hla_threshold = 2
        self.data_test = data_test
        pm = PreprocessedModel(chain)
        dm = Processing(chain)
        df_hla = pm.get_hla_info(hla_threshold,hla_target)
        self.dic_p_train,self.n_pat_train = pm.train_patients_database(hla_target,df_hla)
        self.sign_tcr = pm.get_significant_tcr(hla_target)
        self.df_w_sign_train = pm.look_for_tcr_index(pm.data_train, self.sign_tcr)
        OMs_train, z_train = pm.get_occurrence_matrix_train(self.df_w_sign_train,self.n_pat_train,self.dic_p_train,default_rep)
        self.X_train = pm.get_matches_x(OMs_train,z_train)
        self.y_train = pm.get_matches_y_train(self.n_pat_train,self.dic_p_train,df_hla)
        self.class_weights = mt.get_weights(pm.weights_file,hla_target)
            
        if self.data_test is not None:
            self.dic_p_test,self.n_pat_test = dm.test_patients_database(ps)
            self.df_w_sign_test = pm.look_for_tcr_index(self.data_test, self.sign_tcr)
            if len(self.df_w_sign_test)==0:
                print('No HLA-related TCRs have been found in your repertoire')
                self.hla_prob = 'N/A'
                self.hla_output = False
            else:
                matrix_index = self.renormalize_training_matrix(chain,hla_target,self.df_w_sign_train,self.df_w_sign_test,self.X_train)
                OMs_test, z_test = pm.get_occurrence_matrix_test(OMs_train,self.n_pat_test,self.dic_p_test,matrix_index,self.df_w_sign_test,grouping)
                self.X_test = pm.get_matches_x(OMs_test, z_test)
                if untyped is True:
                    self.hla_prob, self.hla_output = self.get_hla_probabilities(self.X_train,np.array(self.y_train),self.X_test,hla_target,self.class_weights)

    def renormalize_training_matrix(self,chain, hla_target, df_w_sign_train, df_w_sign_test,X_train):
        #l = []
        big_matrix_index = pd.DataFrame()
        for c in chain:
            df_w_sign_train_c = df_w_sign_train.loc[df_w_sign_train['chain']==c,:]
            df_w_sign_train_c.reset_index(drop=True,inplace=True)
            sign_tcr_2 = list(df_w_sign_test['cdr3+v_family'])
            matrix_index = mt.applyParallel(sign_tcr_2, df_w_sign_train_c, mt.get_tcr_index)
            big_matrix_index = pd.concat([big_matrix_index,matrix_index])
            # l0 = list(matrix_index.index)
            # l = l+l0
            
        # if 'alpha' in chain and 'beta' in chain:
        #     l.extend([len(X_train[0])-1])
        # X_norm = X_train[:,l]
        return(big_matrix_index)
    
    def get_hla_probabilities(self,X_train, y_train,X_test,hla_target,class_weight,solver='liblinear',penalty='l1',t=0.9005,l1_ratio=None):
        paramGrid = {
            'C': [50,20,10]+ list(np.logspace(-5,0,50)),
            'solver':[solver], 
            'penalty': [penalty],
            'class_weight':[class_weight]
            }
        lr = LogisticRegression()
        model_lr = GridSearchCV(lr, param_grid=paramGrid, cv=5, verbose=True, n_jobs=-1, scoring='f1_macro')  # neg_log_loss
        
        model_lr.fit(X_train, y_train)
        prob_lr = model_lr.predict_proba(X_test)
        predictions_lr = (prob_lr[0][1] > t) #.astype(int)
        return(prob_lr[0][1],predictions_lr)

    def classifier_params(self,df_param,hla_target):
        
        df1 = df_param.loc[df_param['HLA']==hla_target,:]
        auc = float(df1['AUC'].iloc[0])
        a = float(df1['accuracy'].iloc[0])
        p = float(df1['precision'].iloc[0])
        s1 = float(df1['sensitivity'].iloc[0])
        s2 = float(df1['specificity'].iloc[0])
        return([auc,a,p,s1,s2])
    