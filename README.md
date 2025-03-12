# PRS-Analysis-Nextflow-Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Table of Contents

1.  [Introduction](#Introduction)
2.  [Test-Data](#Test-Data)
3.  [Workflow](#Workflow)
4.  [Installation](#Installation)
5.  [Support](#Support)
6.  [Citation](#Citation)

## Introduction

The pipeline begins by retrieving PGS Score files from the PGS Catalogue, utilizing PGS ids corresponding to various traits. Subsequently, these files are processed to compute PGS scores through PLINK, after which the scores are combined using ElasticNet. Following this consolidation, the pipeline proceeds to generate a prediction score and conduct association analysis specifically tailored to lipid traits.

## Test Data

## Workflow
![pipeline](https://github.com/nanjalaruth/Intergrated_PRS_Analysis/blob/main/conf/pipeline_dag.png)

## Installation 

1. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```
2. R
3. PLINK V1.9
5. Liftover
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
