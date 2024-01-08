import numpy as np
import pandas as pd
import itertools
from functools import partial
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def applyParallel(sign_tcr, df_training, func):

    cores = cpu_count()
    list_split = np.array_split(np.array(sign_tcr), cores)
    func_arg = partial(func, df_training)
    with Pool(cpu_count()) as p:
        ret_list = list(tqdm(p.imap(func_arg, list_split)))
    return pd.concat(ret_list)

def get_tcr_index(df_training,sign_tcr):

    i = []
    for tcr in tqdm(sign_tcr):
        i.append(df_training.index[df_training['cdr3+v_family']==tcr].tolist())
    index = list(itertools.chain(*i))
    df_training.reset_index(drop=True,inplace=True)
    df1 = df_training.iloc[index]
    return(df1)

def get_weights(fInput,hla):
    df = pd.read_csv(fInput,delimiter='\t')
    if hla in list(df['HLA']):
        df_h = df[df['HLA']==hla]
        if df_h['wn'].values > df_h['wp'].values:
            w1 = float(round(df_h.wn/df_h.wp,1))
            w0 = 1
            class_weights = dict([(0, w0), (1, w1)])
        else:
            w0 = round(df_h.wp/df_h.wn,1)
            w1 = 1
            class_weights = dict([(0, w0), (1, w1)])
    else:
        class_weights = None
    return(class_weights)

def print_params(hla,prob,p,i,final_class):
    if prob=='N/A':
        p = ['0.5','1','0','0','1']
    d = {'HLA':[hla], 'Probability':[prob], 'Final_classification':[final_class], 'AUC':[p[0]], 'Accuracy':[p[1]], 'Precision':[p[2]], 'Sensitivity':[p[3]], 'Specificity':[p[4]], 'Individual_ID':[i]}
    df = pd.DataFrame(data=d)
    print('{} '.format(hla)+' -> '+ 'Probability: {}'.format(prob) + '\t Final_classification: {}'.format(final_class)+ '\t AUC: {}'.format(p[0]) + '\t Accuracy: {}'.format(p[1]) + '\t Precision: {}'.format(p[2]) +  '\t Sensitivity: {}'.format(p[3]) + '\t Specificity: {}'.format(p[4]))
    return(df)
