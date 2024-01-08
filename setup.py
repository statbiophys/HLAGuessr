from setuptools import setup, find_packages
import sys
sys.path.append("./HLAGuessr")

def readme():
    with open('README.md') as f:
        return f.read()

data_files_to_include = [('', ['README.md', 'LICENSE','MANIFEST.in'])]

setup(name='HLAGuessr',
      version='0.1.6',
      description='HLA typing from alpha and beta T cell repertoires.',
      long_description='HLAGuessr is a python 3.6 software developed to infer HLA haplotypes from repertoire datasets, using alpha, beta or both chain CDR3 amino acid sequences. The inference is done via a linear model that assigns different weights to a list of TCRs that tend to significantly co-occur (according to Fisher exact test) among people with a given HLA phenotypes. Then, by measuring the presence (or absence) of those TCRs and statistically assessing their concordance with phenotypes of interest, we show that different HLA alleles can be predicted with high accuracy solely on the basis of the TCRβ (and, in some cases,the TCRα)repertoire data generated from peripheral blood.This method was first used in Emerson et al (2017) to study the relation between CMV serostatus and HLA phenotypes from public TCRs. HLAGuessr takes as input a list of TCR CDR3 amino acid sequences and V gene families and optionally a list with the HLA alleles which probabilities will be computed. If no list is specified, a default list of available HLA will be used instead. The output is a list of probabilities of HLA matching followed by the parameters (AUC, accuracy, precision, sensitivity and specificity) measured over an external validation dataset that can be used as guidance for the certitude of the computation',
      url='https://github.com/statbiophys/HLAGuessr',
      author='María Ruiz Ortega',
      author_email='ruizormaria@gmail.com',
      license='GPLv3',
      classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Intended Audience :: Healthcare Industry',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Programming Language :: Python :: 3.6',
            ],
      packages=find_packages(),
      install_requires=['numpy','scikit-learn','pandas','tqdm'], # 'multiprocessing',
      package_data = {
            'Training_data': ['Training_data/*'],
            'Validation_data': ['Validation_data/*.tsv'],
            },
      data_files = data_files_to_include,
      include_package_data=True,
      entry_points = {'console_scripts': [
            'HLAGuessr-infer_hla=HLAGuessr.infer_hla:main',], },
      zip_safe=False)