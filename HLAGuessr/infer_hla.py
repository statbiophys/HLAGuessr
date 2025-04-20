from __future__ import print_function, division,absolute_import
from HLAGuessr.processing import Processing   
from HLAGuessr.evaluate_data import GetProbabilities
import HLAGuessr.model_tools as mt
import pandas as pd, numpy as np
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from optparse import OptionParser
        
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line script to evaluate the probability of HLA allele presence.

    Copyright (C) 2023 María Ruiz Ortega

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

def main():
    """ Evaluate repertoires."""
    parser = OptionParser(conflict_handler="resolve")

    # input output
    parser.add_option('-a', '--alpha_infile', dest = 'alpha_infile',metavar='PATH/TO/FILE', help='read in CDR3 alpha sequences and V gene families from PATH/TO/FILE')
    parser.add_option('-b', '--beta_infile', dest = 'beta_infile',metavar='PATH/TO/FILE', help='read in CDR3 alpha sequences and V gene families from PATH/TO/FILE')
    parser.add_option('-o', '--outfile', dest = 'outfile_name', metavar='PATH/TO/FILE', help='write HLA probabilities to PATH/TO/FILE')
    
    # delimeters
    parser.add_option('-d', '--delimiter', type='choice',default=None, dest='delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare infile delimiter. Default is tab for .tsv input files, comma for .csv files, and any whitespace for all others. Choices: 'tab', 'space', ',', ';', ':'")
    
    # HLA alleles
    parser.add_option('--hla', dest='hla', default=False, metavar=' "X*aa:bb,Y*cc:dd,..." ', help="specify with quotation marks the HLA allele (or HLA alleles, separated by ,) whose probabilities want to be inferred. If nothing is specified, a default list of alleles will be used.")
    
    # parser.add_option('--seq_in', '--seq_index', type='int', metavar='INDEX', dest='seq_in_index', default = 0, help='specifies sequences to be read in are in column INDEX. Default is index 0 (the first column).')
    # parser.add_option('-m', '--max_number_of_seqs', type='int',metavar='N', dest='max_number_of_seqs', help='evaluate for at most N sequences.')
    # parser.add_option('--lines_to_skip', type='int',metavar='N', dest='lines_to_skip', default = 0, help='skip the first N lines of the file. Default is 0.')

    (options, args) = parser.parse_args()

    # Load repertoire files
    
    chain = []

    if options.alpha_infile is not None:
        infile_name = options.alpha_infile
        chain.append('alpha')
        
        if not os.path.isfile(options.alpha_infile):
            print('Cannot find alpha input file: ' + infile_name)
            print('Exiting...')
            return(-1)
        
    if options.beta_infile is not None:
        infile_name = options.beta_infile
        chain.append('beta')
        
        if not os.path.isfile(options.beta_infile):
            print('Cannot find beta input file: ' + infile_name)
            print('Exiting...')
            return(-1)

    # Parse delimiter
    delimiter = options.delimiter
    if delimiter is None: #Default case
        if options.infile_name is None:
            delimiter = '\t'
        elif infile_name.endswith('.tsv'): # parse TAB separated value file
            delimiter = '\t'
        elif infile_name.endswith('.csv'): # parse COMMA separated value file
            delimiter = ','
    else:
        try:
            delimiter = {'tab': '\t', 'space': ' ', ',': ',', ';': ';', ':': ':'}[delimiter]
        except KeyError:
            pass # Other string passed as the delimiter.
        
    
    process = Processing(chain,options.alpha_infile,options.beta_infile,delimiter)    
    main_path = os.getcwd()
    hla_threshold=1
    grouping = False
    untyped = True
    f_params = main_path + '/Training_data/classifier_params_alpha+beta.tsv'
    df_param = pd.read_csv(f_params,delimiter='\t')
    
    hla_set = []
    
    print(options.hla)
    
    if options.hla is not False:
        hla_0 = options.hla.split(',')
        for h in hla_0:
            if h not in list(df_param['HLA']):
                print('No available information for HLA '+h)
            else:
                hla_set.append(h)

    else:
        hla_set = list(df_param['HLA'])
    
    df_final_params = pd.DataFrame()
    
    if len(hla_set)==0:
            print('No available information for the HLA alleles provided')
            print('Exiting...')
            return(-1)
         
    else:
        
        non_empty = [p for p in process.ps if p]
        unique_ps = list({tuple(p) for p in non_empty})
        unique_ps = [list(p) for p in unique_ps]

        for p in unique_ps:
            
            print('** Individual ' +p[0]+ ': **')

            for hla_target in hla_set:
                
                ev = GetProbabilities(hla_target,chain,hla_threshold,
                                      process.ps,process.data_test,grouping,untyped)
                params = ev.classifier_params(df_param,hla_target)
                
                df = mt.print_params(hla_target,ev.hla_prob,params,p[0],ev.hla_output)
                df_final_params = pd.concat([df,df_final_params],ignore_index=True)

            if options.outfile_name is not None:
                outfile_name = options.outfile_name
                df_final_params.to_csv(outfile_name,sep='\t',index=False)
                print('Output '+ outfile_name + ' succesfully generated.')

if __name__ == '__main__': main()