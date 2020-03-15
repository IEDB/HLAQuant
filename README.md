# HLAQuant
### Author: Austin Crinklaw

## What is it?
HLAQuant is a pipeline that produces fast and accurate allele specific expression for HLA genes. This is done by using sequences specific to a donor's HLA typing to quantify expression.

## Requirements:
- Linux OS
- NCBI Blast+
- Salmon -- please ensure Salmon is on your PATH environment variable.
- Python 3+
  - Python packages: Pandas, BioPython

## How to use:

### Installation:
HLAQuant can be downloaded through PyPI using the following pip command.
```shell
pip install hlaquant
```

### Input  
HLAQuant takes two input files currently. 
- File one (-hla) consists of a sample_id and then a list of alleles corresponding to that sample's typing (tab separated)
- File two (-fastq) consists of a sample_id and then a list of FASTQ files corresponding to that sample (paired-end or single-end)
Examples of these inputs can be found under the 'test_data/' directory

### Usage
- A list of parameters and their descriptions can be found with the -h flag
```shell
python -m HLAQuant -h
```

### Output  
The output will match that of Salmons. It consists of a tab separated file containing the transcript ID (in this case, a specific HLA allele), as well as the number of reads. TPMs can be ignored as they will be inaccurate since we are only quantifying over a few sequences.

## How does it work?
- We first take the list of alleles and fetch the corresponding sequences from IMGT
- We build an index for quantification using these sequences
- We then perform quantification using this personalized index


## References:
This pipeline would be unable to work without Salmon.

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.

More recently, Aguiar et al. display evidence of improving expression estimates and eQTL analysis using personalized sequences

Aguiar VRC, César J, Delaneau O, Dermitzakis ET, Meyer D. Expression estimation and eQTL mapping for HLA genes with a personalized pipeline. PLoS Genet. 2019;15(4):e1008091. Published 2019 Apr 22. doi:10.1371/journal.pgen.1008091

