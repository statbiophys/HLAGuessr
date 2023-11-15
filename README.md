# HLAGuessr
HLAGuessr is a python 3.6 software developed to infer HLA haplotypes from repertoire datasets, using alpha, beta or both chain CDR3 amino acid sequences. The inference is done via a linear model that assigns different weights to a list of TCRs that tend to significantly co-occur (according to Fisher exact test) among people with a given HLA phenotypes. Then, by measuring the presence (or absence) of those TCRs and statistically assessing their concordance with phenotypes of interest, we show that different HLA alleles can be predicted with high accuracy solely on the basis of the TCRβ (and, in some cases,the TCRα)repertoire data generated from peripheral blood.This method was first used in Emerson et al. (2017) to study the relation between CMV serostatus and HLA phenotypes from public TCRs. HLAGuessr takes as input a list of TCR CDR3 amino acid sequences and V gene families and optionally a list with the HLA alleles which probabilities will be computed. If no list is specified, a default list of available HLA will be used instead. The output is a list of probabilities of HLA matching followed by the parameters (AUC, accuracy, precision, sensitivity and specificity) measured over an external validation dataset that can be used as guidance for the certitude of the computation

## Contact

Any issues or questions should be addressed to [us](mailto:ruizormaria@gmail.com).

## License

Free use of HLAGuessr is granted under the terms of the GNU General Public License version 3 (GPLv3).

## References 

[1] Emerson RO, et al. (2017) Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA- mediated effects on the T cell repertoire. Nature Genetics 49:659–665.
