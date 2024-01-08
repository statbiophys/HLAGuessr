# HLAGuessr
HLAGuessr is a python 3.6 software developed to infer HLA haplotypes from repertoire datasets, using alpha, beta or both chain CDR3 amino acid sequences. The inference is done via a linear model that assigns different weights to a list of TCRs that tend to significantly co-occur (according to Fisher exact test) among people with a given HLA phenotypes. Then, by measuring the presence (or absence) of those TCRs and statistically assessing their concordance with phenotypes of interest, we show that different HLA alleles can be predicted with high accuracy solely on the basis of the TCRβ (and, in some cases,the TCRα)repertoire data generated from peripheral blood.This method was first used in Emerson et al. (2017) to study the relation between CMV serostatus and HLA phenotypes from public TCRs. HLAGuessr takes as input a list of TCR CDR3 amino acid sequences and V gene families and optionally a list with the HLA alleles which probabilities will be computed. If no list is specified, a default list of available HLA will be used instead. The output is a list of probabilities of HLA matching followed by the parameters (AUC, accuracy, precision, sensitivity and specificity) measured over an external validation dataset that can be used as guidance for the certitude of the computation.

## Version
Latest released version: 0.1.6


## Installation
HLAGuessr is available on PyPI and can be downloaded and installed through pip:


 ```pip install HLAGuessr```

HLAGuessr is also available on [GitHub](https://github.com/statbiophys/HLAGuessr). The command line entry points can be installed by using the setup.py script:

 ```$python setup.py install .```.

 Directory architecture:
```
 HLAGuessr/
│   README.md
│   LICENSE
│   setup.py
│   MANIFEST.in  
│
└─── HLAGuessr/
    │   __init__.py
    │   processing.py
    │   evaluate_data.py
    │   infer_hla.py
    │   load_model.py
    │   model_tools.py
    │
    └─── Training data/
        │    alpha_cdr3aa_shared_at_least_3.txt 
        │    beta_cdr3aa_shared_at_least_3.txt
        │    classifier_params_alpha+beta.tsv
        │    classifier_params_alpha.tsv
        │    classifier_params_beta.tsv
        │    HLA_grouped_patients.tsv
        │    inference_all_data_TCRs_sign_HLA.tsv
        │
        
    └─── Example_validations_data/
        │    alpha_example.tsv 
        │    beta_example.tsv
```

## Command line console script

There is one command line console scripts (the scripts can still be called as executables if HLAGuessr is not installed):
1. HLAGuessr-infer_hla
  * It computes the probability that a list of HLA alleles is present in alpha and/or repertoire data

You can execute it with the -h or --help flags to get the options.

### General commands summary

| Selected Options                               | Description                                      |
|------------------------------------------------|--------------------------------------------------|
|   **-h**, **--help**                           |   show full Options list and exit                          |
|   **-a**, **--alpha_infile** PATH/TO/ALPHA_FILE|   read input alpha repertoire file containing CDR3α sequences from PATH/TO/ALPHA_FILE |
|   **-b**, **--beta_infile** PATH/TO/BETA_FILE  |   read input beta repertoire file containing CDR3β sequences from PATH/TO/BETA_FILE |
|   **-o**, **--outfile** PATH/TO/FILE           |   write probabilities and report of the classifier parameters in PATH/TO/FILE |
|   **--hla** HLA ALLELE (or list of alleles)    |   specify HLA allele (or list of them separated for a comma) to compute probabilities. If not allele is specified, a default list of available alleles will be entirely explored.  |
|  **-d** DELIMITER                              |   declare infile delimiter. Default is tab for .tsv input files, comma for .csv files, and any whitespace for all others. Choices: 'tab', 'space', ',', ';', ':'    |

### Notes about input files format

The script has only minimal file parsing built in, so reading in sequences from a file requires the file to be structured in a particular way. Data must be presented with delimiter spaced values (i.e. the data is organized in columns separated by delimiter like a tab for .tsv or a comma for .csv file). The first column, called *cdr3aa*, must contain the **CDR3 region** in amino acid format, from the conserved cysteine C (INCLUSIVE) in the V region to the conserved F (INCLUSIVE) in the J; the second column, called *v_family*, must contain the **V gene family** in the format TR(A/B)VX; the third column, named *Patient*, must contain some **string identifier for each individual** included in the repertoire data. An example of how it should be presented:
```
    cdr3aa         v_family	   Patient
CAAAADAGGTSYGKLTF   TRAV23  	    B1
CAAAAFGNEKLTF       TRAV29  	    B1
CAAAAGANNLFF        TRAV23  	    B1
CAAAAGGTSYGKLTF     TRAV29  	    B1
```
## Output

As a result for the computation, for each patient and HLA allele, the classifier will print a probability of the allele to be present in the person and a summary of the parameters that caracterized the classifier performance (AUC, accuracy, precision sensitivity and specificity) when it was applied to an external validation dataset. By displying these metrics, the user can weigh the validity of the prediction. In order to reduce the number of false positives, only those alleles with probability higher than 0.9 will be consider positives for the analysed repertoire.
If a file is specified to write to (using -o, see Options), the generated predictions are written to the file, otherwise they are printed to stdout.

## Quick Demo
After installing HLAGuessr, we offer a quick demonstration of the console scripts using two example files that can be found in Example_validation_data folder, alpha_example and beta_example. 
```
1. HLAGuessr-infer_hla -a ~/HLAGuessr/HLAGuessr/Example_validation_data/alpha_example.tsv -b ~/HLAGuessr/HLAGuessr/Example_validation_data/beta_example.tsv --hla "A*02:01" -d tab

100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 1842.48it/s]
100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 2011.12it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 2195.40it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1759.91it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1991.48it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1763.24it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1780.93it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1955.84it/s]
8it [00:00, 330.90it/s]
A*02:01  -> Probability: 1.0	 AUC: 0.9879227053140096	 Accuracy: 0.7073170731707317	 Precision: 0.6	 Sensitivity: 1.0	 Specificity: 0.4782608695652174

2. HLAGuessr-infer_hla -a ~/HLAGuessr/HLAGuessr/Example_validation_data/alpha_example.tsv -b ~/HLAGuessr/HLAGuessr/Example_validation_data/beta_example.tsv --hla "B*02:01" -d tab
No available information for HLA B*02:01
No available information for the HLA alleles provided
Exiting...

3. HLAGuessr-infer_hla -a ~/Scripts/HLA_Guessr_package/HLAGuessr/Example_validation_data/alpha_example.tsv -b ~/Scripts/HLA_Guessr_package/HLAGuessr/Example_validation_data/beta_example.tsv --hla "A*02:01,B*07:02" -d tab
Individual 17_B:
100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 2141.29it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 2096.76it/s]
100%|███████████████████████████████████████████| 9/9 [00:00<00:00, 1864.32it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1885.72it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1767.05it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1716.78it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1910.52it/s]
100%|███████████████████████████████████████████| 8/8 [00:00<00:00, 1881.49it/s]
8it [00:00, 357.27it/s]
A*02:01  -> Probability: 1.0	 AUC: 0.9879227053140096	 Accuracy: 0.7073170731707317	 Precision: 0.6	 Sensitivity: 1.0	 Specificity: 0.4782608695652174
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2317.84it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2303.79it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2086.16it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2076.94it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2227.29it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2659.19it/s]
100%|█████████████████████████████████████████| 28/28 [00:00<00:00, 2137.85it/s]
100%|█████████████████████████████████████████| 27/27 [00:00<00:00, 2323.19it/s]
8it [00:00, 213.73it/s]
B*07:02  -> Probability: 1.0	 AUC: 1.0	 Accuracy: 0.7804878048780488	 Precision: 0.5263157894736842	 Sensitivity: 1.0	 Specificity: 0.7096774193548387
```


## Contact

Any issues or questions should be addressed to [us](mailto:ruizormaria@gmail.com).

## License

Free use of HLAGuessr is granted under the terms of the GNU General Public License version 3 (GPLv3).

## References 

[1] Emerson RO, et al. (2017) Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA- mediated effects on the T cell repertoire. Nature Genetics 49:659–665.
