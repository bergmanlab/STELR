<p align="center">
    <img src="https://github.com/bergmanlab/TELR/blob/master/TELR.png?raw=true" alt="TELR"/>
</p>

## TELR
TELR (pronounced Teller) is a fast non-reference transposable element (TE) detector from long read sequencing data (PacBio or Oxford Nanopore).

## Installation
You can use conda to install dependencies and create running environment for TELR.
1. install all dependencies
```
conda create -n TELR_env
conda activate TELR_env
conda install -y biopython
conda install -c conda-forge -y pandas
conda install -y repeatmasker=4.0.7
conda install -y samtools=1.9
conda install -y bcftools=1.9
conda install -y bedtools
conda install -y ngmlr=0.2.7
conda install -y sniffles=1.0.11
conda install -y wtdbg
conda install -y seqtk
conda install -y minimap2

```
2. Install TELR
```
git clone git@github.com:bergmanlab/TELR.git
cd TELR
python3 main.py --help</pre>
```

## Parameters
```
usage: main.py [-h] -i READ -r REFERENCE -l LIBRARY [-x PRESETS] [-o OUT]
               [-t THREAD] [-g GAP] [-p OVERLAP]

Script to detect TEs in long read data

required arguments:
  -i READ, --read READ
      reads in fasta/fastq format
  -r REFERENCE, --reference REFERENCE
      reference genome in fasta format
  -l LIBRARY, --library LIBRARY
      TE consensus sequences in fasta format

optional arguments:
  -h, --help
      show this help message and exit
  -x PRESETS, --presets PRESETS
      parameter presets for different sequencing technologies
  -o OUT, --out OUT
      directory to output data (default = '.')
  -t THREAD, --thread THREAD
      max cpu threads to use (default = '1')
  -g GAP, --gap GAP
      max gap size for flanking sequence alignment (default = '20')
  -p OVERLAP, --overlap OVERLAP
      max overlap size for flanking sequence alignment (default = '20')
```