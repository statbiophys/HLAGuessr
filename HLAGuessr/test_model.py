#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd, numpy as np, HLAGuessr.model_tools as mt
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import (roc_curve,ConfusionMatrixDisplay,confusion_matrix)
from HLAGuessr.load_model import PreprocessedModel
from HLAGuessr.evaluate_data import GetProbabilities


class CrossValidation(object):
    
    def __init__(self,hla_target,hla_input,chain,ps,data_test,grouping,hla_threshold,solver=None,penalty=None,untyped=False):
        
        if hla_threshold == None:
            hla_threshold = 2
        
        ev = GetProbabilities(hla_target,chain,hla_threshold,ps,data_test,grouping,untyped)
        pm = PreprocessedModel(chain,hla_file=hla_input)
        df_hla = pm.get_hla_info(hla_threshold,hla_target)
        self.y_test = self.get_matches_y_test(ev.n_pat_test,ev.dic_p_test,df_hla)
        class_weights = mt.get_weights(pm.weights_file,hla_target)
        if solver is None:
            solver='liblinear'
        if penalty is None:
            penalty = 'l1'
        self.df_summary, self.model_lr = self.learn_and_plot(ev.X_train,ev.y_train,ev.X_test,self.y_test,chain,hla_target,class_weights,solver,penalty)
        
    def get_matches_y_test(self,n_patients,dic,df_hla):
        
        OH = np.zeros([1,int(n_patients)])
        for i,x in enumerate(df_hla['pos_patients'].values):
            for c in str(x).split(', '):
                if c in dic:
                    i = dic.get(str(c)) 
                    OH[0][int(i)]=1
        y = list(OH[0])
        return(y)
    
    def learn_and_plot(self,X_train,y_train,X_test,y_test,chain,hla_target,class_weights,solver,penalty,t=0.9005,l1_ratio=None,plot=False,directory=None):
        
        auc_t, accuracy, precision, recall, specificity = [],[],[],[],[]
        paramGrid = {
            'C': [50,20,10]+ list(np.logspace(-5,0,50)),
            'solver':[solver], 
            'penalty': [penalty],
            'class_weight':[class_weights]
            }
        lr = LogisticRegression()
        model_lr = GridSearchCV(lr, param_grid = paramGrid, cv = 5, verbose=True, n_jobs=-1, scoring='f1')
        model_lr.fit(X_train, y_train)
        prob_lr = model_lr.predict_proba(X_test)[:,1]
        predictions_lr = (prob_lr > t).astype(int)

        cnf_matrix = confusion_matrix(y_test, predictions_lr)
        auc_t.append(metrics.roc_auc_score(y_test, prob_lr))
        accuracy.append(metrics.accuracy_score(y_test, predictions_lr))
        precision.append(metrics.precision_score(y_test, predictions_lr))
        recall.append(metrics.recall_score(y_test, predictions_lr))
        specificity.append(cnf_matrix[0][0]/(cnf_matrix[0][1]+cnf_matrix[0][0]))
        fpr, tpr, thresholds = roc_curve(y_test, model_lr.predict_proba(X_test)[:,1])

        #### Model summary ####
        #print(classification_report(y_test, predictions_lr))

        #### Plot model summary ####
        if plot==True:
            display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=metrics.roc_auc_score(y_test, prob_lr))
            display.plot()
            if directory!=None:
                plt.savefig(directory+'roc_curve.png')

            ConfusionMatrixDisplay.from_predictions(y_test, predictions_lr)
            if directory!=None:
                plt.savefig(directory+'confusion_matrix.png')

        cols = ['HLA','AUC','accuracy','precision','sensitivity','specificity']
        df = pd.DataFrame(list(zip([hla_target],auc_t,accuracy,precision,recall,specificity)),columns=cols)
        if len(chain)==2:
            df['chain'] = 'alpha+beta'
        else:
            df['chain'] = chain[0]

        return(df, model_lr)
    