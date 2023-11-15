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
    df1 = df_training.iloc[index]
    return(df1)

def print_params(hla,prob,p):
    if prob=='N/A':
        p = ['N/A','N/A','N/A','N/A','N/A']
    d = {'HLA': [hla], 'Probability': [prob], 'AUC': [p[0]], 'Accuracy':[p[1]], 'Precision':[p[2]], 'Sensitivity':[p[3]], 'Specificity':[p[4]]}
    df = pd.DataFrame(data=d)
    print('{} '.format(hla)+' -> '+ 'Probability: {}'.format(prob) + '\t AUC: {}'.format(p[0]) + '\t Accuracy: {}'.format(p[1]) + '\t Precision: {}'.format(p[2]) +  '\t Sensitivity: {}'.format(p[3]) + '\t Specificity: {}'.format(p[4]))
    return(df)
