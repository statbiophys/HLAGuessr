# HLAGuessr
HLAGuessr is a python 3.6 software developed to infer HLA haplotypes from repertoire datasets, using alpha, beta or both chain CDR3 amino acid sequences. The inference is done via a linear model that assigns different weights to a list of TCRs that tend to significantly co-occur (according to Fisher exact test) among people with a given HLA phenotypes. Then, by measuring the presence (or absence) of those TCRs and statistically assessing their concordance with phenotypes of interest, we show that different HLA alleles can be predicted with high accuracy solely on the basis of the TCRβ (and, in some cases,the TCRα)repertoire data generated from peripheral blood.This method was first used in Emerson et al. (2017) to study the relation between CMV serostatus and HLA phenotypes from public TCRs. HLAGuessr takes as input a list of TCR CDR3 amino acid sequences and V gene families and optionally a list with the HLA alleles which probabilities will be computed. If no list is specified, a default list of available HLA will be used instead. The output is a list of probabilities of HLA matching followed by the parameters (AUC, accuracy, precision, sensitivity and specificity) measured over an external validation dataset that can be used as guidance for the certitude of the computation.

## Version
Latest released version: 0.0.2

## Installation
HLAGuessr is available on PyPI and can be downloaded and installed through pip:

 ```pip install HLAGuessr==0.0.2```.

HLAGuessr is also available on [GitHub](https://github.com/mariaruizortega94/HLAGuessr). The command line entry points can be installed by using the setup.py script:

 ```$python setup.py install .```.

 Directory architecture:

 HLAGuessr/
│   README.md
│   LICENSE
│   setup.py
│   MANIFEST.in  
│
└───HLAGuessr/
    │   __init__.py
    │   processing.py
    │   evaluate_data.py
    │   infer_hla.py
    │   load_model.py
    │   model_tools.py
    │
    └───Training data/
        │    alpha_cdr3aa_shared_at_least_3.txt 
        │    beta_cdr3aa_shared_at_least_3.txt
        │    classifier_params_alpha+beta.tsv
        │    classifier_params_alpha.tsv
        │    classifier_params_beta.tsv
        │    HLA_grouped_patients.tsv
        │    inference_all_data_TCRs_sign_HLA.tsv
        │
        
    └───Example_validations_data/
        │    alpha_example.tsv 
        │    beta_example.tsv


## Command line console scripts and Examples

There is one command line console scripts (the scripts can still be called as executables if HLAGuessr is not installed):
1. HLAGuessr-infer_hla
  * It computes the probability that a list of HLA alleles is present in alpha and/or repertoire data

You can execute it with the -h or --help flags to get the options.

### Quick Demo
After installing HLAGuessr, we offer a quick demonstration of the console scripts using two example files that can be found in Example_validation_data folder, alpha_example and beta_example. 
1. ```HLAGuessr-infer_hla -a ~/HLAGuessr/HLAGuessr/Example_validation_data/alpha_example.tsv -b ~/HLAGuessr/HLAGuessr/Example_validation_data/beta_example.tsv --hla A*02:01 -d tab

100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 1842.48it/s]
100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 2011.12it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 2195.40it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1759.91it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1991.48it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1763.24it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1780.93it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1955.84it/s]
8it [00:00, 330.90it/s]
A*02:01  -> Probability: 1.0	 AUC: 0.9879227053140096	 Accuracy: 0.7073170731707317	 Precision: 0.6	 Sensitivity: 1.0	 Specificity: 0.4782608695652174```


## Contact

Any issues or questions should be addressed to [us](mailto:ruizormaria@gmail.com).

## License

Free use of HLAGuessr is granted under the terms of the GNU General Public License version 3 (GPLv3).

## References 

[1] Emerson RO, et al. (2017) Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA- mediated effects on the T cell repertoire. Nature Genetics 49:659–665.
