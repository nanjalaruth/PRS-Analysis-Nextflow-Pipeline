# PRS-Analysis-Nextflow-Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Table of Contents

1.  [Introduction](#Introduction)
2.  [Workflow](#Workflow)
3.  [Installation](#Installation)
4.  [Support](#Support)
5.  [Citation](#Citation)

## Introduction

The pipeline begins by retrieving PGS Score files from the PGS Catalogue, utilizing PGS ids corresponding to various traits. Subsequently, these files are processed to compute PGS scores through PLINK, after which the scores are combined for each trait in preparation for further analysis using ElasticNet

## Workflow
![pipeline](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/output.png)

## Installation 
### Data
1. SNP Info files
The info files should be in a tsv format with header information, as shown in the examples below. If your file is in any other format, please convert it to this format using bash, python or R
The examples are:
a. hg38 info files named as hg38_snp_info_header.txt

| chr_name  | chr_position | rsID    |
|-----------|--------------|---------|
| 10   | 76684698    | rs241         |
| 10   | 96480625    | rs243         |
| 10   | 20703742    | rs244         |  

b. hg37 info files named as hg37_snp_info_header.txt
| chr_name  | chr_position | rsID    |
|-----------|--------------|---------|
| 10   | 76684698    | rs241         |
| 10   | 96480625    | rs243         |
| 10   | 20703742    | rs244         |  

2. The liftover chain files
Please download:
a. hg19ToHg38.over.chain.gz
b. hg38ToHg19.over.chain.gz

3. Have your gentype data ready in PLINK format

### Tools
1.Conda

a. [Download Miniconda](https://www.anaconda.com/download/) for your specific OS to your home directory
    - Linux: `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
    - Mac: `curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`
b. Run:
    - `bash Miniconda3-latest-Linux-x86_64.sh`
    - `bash Miniconda3-latest-MacOSX-x86_64.sh`
c. Follow all the prompts: if unsure, accept defaults
d. Close and re-open your terminal
e. If the installation is successful, you should see a list of installed packages with
    - `conda list`
If the command cannot be found, you can add Anaconda bin to the path using:
    ` export PATH=~/miniconda3/bin:$PATH`
    
2. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```

3. R
```
conda install R
```
   
5. PLINK V1.9
```
https://www.cog-genomics.org/plink/
```
   
6. Liftover
We used conda for installation but feel free to install the liftover tool using other options
```
#Installation using conda
conda create -n liftover
conda activate liftover
conda install -c bioconda ucsc-liftover
```

## Running the pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

### Own data
Start running your own analysis either by using flags as shown below:

Run pipeline directly by providing paths to the input files

Alternatively
 Run your own analysis by modifying the conf/test.config file to suit the path to your data location and then run the command as below:
 
 ```
 nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline -profile slurm -c <path to your edited config file> -resume
 ```

## To run the updated version of this pipeline, run:

 ```
 nextflow pull nanjalaruth/Intergrated_PRS_Analysis
 ```
## Output
[Analysis_output](https://nanje.quarto.pub/intergrated_prs/)
[Nextflow report]()

## Support
I track open tasks using github's [issues](https://github.com/nanjalaruth/Intergrated_PRS_Analysis/issues)

## Citation
