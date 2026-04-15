# PRS-Analysis-Nextflow-Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)


## Table of Contents
1.  [Introduction](#Introduction)
2.  [Installation](#Installation)
3.  [Running the pipeline](#Running-the-pipeline)
4.  [Workflow](#Workflow)
5.  [Output](#Output)
6.  [Support](#Support)
7.  [Citation](#Citation)


## Introduction
In the era of large-scale genomics, efficiently computing Polygenic Scores (PGS) across multiple phenotypes and score IDs is a critical yet complex task. Manual processing is not only time-consuming but also prone to errors, making it difficult to ensure reproducibility and scalability. Our `Nextflow` pipeline automates the entire PGS computation workflow, enabling seamless integration of genotype data, PGS weights, and phenotype information. By leveraging parallelization, error handling, and robust quality control, this pipeline ensures that PGS scores are computed accurately and efficiently across diverse datasets. Designed for scalability, it allows researchers to process multiple phenotypes with multiple PGS ids simultaneously, making large-scale genetic studies more reproducible, efficient, and easy to maintain.

The pipeline begins by retrieving PGS Score files from the PGS Catalogue, utilizing PGS IDs corresponding to various traits. Before computing PGS scores, the files are processed to generate necessary inputs, including modifications for liftover (genomic coordinate conversion) and formatting for PLINK-based PGS score calculation. Once the scores are computed, they are combined for each trait in preparation for further analysis using ElasticNet.


## Installation 
### Tools

`The pipeline requires Nextflow, PLINK, and R to run`
To install the tools follow the commands below:

1. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```

2. PLINK V1.9
```
https://www.cog-genomics.org/plink/
```

3. Conda

`conda will help with installation of R. If you already have R installed, you can skip Miniconda installation`

- [Download Miniconda](https://www.anaconda.com/download/) for your specific OS to your home directory
    - Linux: `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`
    - Mac: `curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`
- Run:
    - `bash Miniconda3-latest-Linux-x86_64.sh`
    - `bash Miniconda3-latest-MacOSX-x86_64.sh`
- Follow all the prompts: if unsure, accept defaults
- Close and re-open your terminal
- If the installation is successful, you should see a list of installed packages with
    - `conda list`
- If the command cannot be found, you can add Anaconda bin to the path using:
    `export PATH=~/miniconda3/bin:$PATH`
    
4. R
```
conda install R
```

### Data
1. PGS IDs files
- Create a folder on your local directory called "all_blood_traits_prs_scores", which should contain all PGS ids files with suffix `*_PGS_score_ids.txt`. The prefix for the files should be similar to the `bloodCells` in the step `#Running-the-pipeline`.

    - For example:
        - If you are interested in the phenotype, `basophil`, the folder should have a PGS score file called `baso_PGS_score_ids.txt`. 
        - The content of the file `baso_PGS_score_ids.txt` should ge PGS ids for basophils copied from [PGS catalogue](https://www.pgscatalog.org/) in the format shown below.  
          ||
          |-----------|
          |PGS003940|
          |PGS004727|
          |PGS004728|
      
        - PS: The file should have no header, just the IDS.
 
 `Note: The pipeline will download data in the GRCh37/B37 format only`

2. Have your genotype data ready in PLINK format ie bed, bim, fam


##  Running-the-pipeline
### Required Arguments

| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| --basePath  | /path/to/all_blood_traits_prs_scores | Directory pattern for PGS Ids files      | 
| --bloodCells  | baso\rbc\wbc | Names of your phenotypes | 
| --plink_file  | bed\bim\fam | PLINK genotype files   | 

- The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.
- Modify the conf/test.config file particularly the lines below to suit the path to your data location:
```
params.bloodCells = ["baso", "rbc", "wbc"]
params.basePath = "/path/to/all_blood_traits_prs_scores/folder"

plink_file = [
    ['UGRC', '/path/to/uganda.bed', '/path/to/uganda.bim', '/path/to/uganda.fam']
]
```

Then run command as below:
```
nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline -resume -c <path to your edited conf/test.config file> 
```

If you have cloned the Github repository, and modified the config file on your terminal, run it as:
```
nextflow run main.nf -resume -c conf/test.config
```

### To run the updated version of this pipeline, run:

 ```
 nextflow pull nanjalaruth/PRS-Analysis-Nextflow-Pipeline
 ```

## Workflow
A summary of the steps followed in our analysis include;
- Downloading score files from PGS catalogue
    - modify_score_file (Removes header from output file)
- Downloading metadata files from PGS catalogue
    - modify metadata file format   
- modify_score_file_3 (makes sure there is a SNP identifier, and outputs a clean three-column score file suitable for polygenic scoring)
- Computes PGS scores using `PLINK V1.9`
- Modify the output and concatenate scores for each phenotype
  
![pipeline](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/output.png)

## Output
[Analysis_output](https://nanje.quarto.pub/intergrated_prs/)
[Nextflow report](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/execution_report.html)

## Support
I track open tasks using github's [issues]([https://github.com/nanjalaruth/IPRS-Analysis-Nextflow-Pipeline/issues](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/issues))

## Citation
