# HLAQuant
### Author: Austin Crinklaw

## What is it?
HLAQuant is a pipeline that produces fast and accurate allele specific expression for HLA genes. This is done by quantifying on the peptide binding groove domain using Salmon with personalized sequences.

## Requirements:
- Linux OS
- NCBI Blast+
- Salmon (distributed with package)
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
- Next we extract the sequences corresponding to their groove domains from these sequences
- We build an index for quantification using these G-domain sequences
- We then perform quantification using this index

The paper outlining this method in detail can be found [....somewhere when it is published]

## References:
This pipeline would be unable to work without Salmon
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
