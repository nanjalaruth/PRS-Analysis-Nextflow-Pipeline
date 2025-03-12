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
### Data
1. PGS IDs and Phenotype files
- Create a folder on your local directory called "all_blood_traits_prs_scores", which should contain all PGS ids files with suffix `*_PGS_score_ids.txt` and phenotype files with the suffix `*_pheno.tsv`. The prefix for the files should be the same ones used as a tag for the `bloodCells` as used in the subsequent step ` [Running the pipeline](#Running-the-pipeline)` e.g. `baso`.
- For example if you are interested in the phenotype, `basophil`, the folder should have:
  - Phenotype file named `baso_pheno.tsv` in this format

      | FID  | IID | baso    |
      |-----------|--------------|---------|
      |94343      |94343         |1.1|
      |94326      |94326         |0.9|
      |94097      |94097         |0.8|

  - PGS score ids called `baso_PGS_score_ids.txt`. 
    - Copy and paste ids from PGS catalogue in the format shown below.
    - The file should have no header, just the IDS.
  
      ||
      |-----------|
      |PGS003940|
      |PGS004727|
      |PGS004728|

2. SNP Info files for hg19/GRCh37 and hg38/GRCh38 genome builds
- The info files should be in a tsv format with header information, as shown in the examples below.
- If your file is in any other format, please convert it to this format using bash, python or R
- The examples are:
  - hg38 info files named as hg38_snp_info_header.txt

    | chr_name  | chr_position | rsID    |
    |-----------|--------------|---------|
    | 10   | 76684698    | rs241         |
    | 10   | 96480625    | rs243         |
    | 10   | 20703742    | rs244         |  

  - hg37 info files named as hg37_snp_info_header.txt
    | chr_name  | chr_position | rsID    |
    |-----------|--------------|---------|
    | 10   | 78444456    | rs241         |
    | 10   | 98240382     | rs243         |
    | 10   | 20992671    | rs244         |  

3. The liftover chain files
- Please download the chain files:
  - hg19ToHg38.over.chain.gz
  - hg38ToHg19.over.chain.gz

4. Have your gentype data ready in PLINK format

### Tools
1.Conda
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
    ` export PATH=~/miniconda3/bin:$PATH`
    
2. Nextflow
```
wget -qO- https://get.nextflow.io | bash
```

3. R
```
conda install R
```
   
4. PLINK V1.9
```
https://www.cog-genomics.org/plink/
```
   
5. Liftover
We used conda for installation but feel free to install the liftover tool using other options
```
conda create -n liftover
conda activate liftover
conda install -c bioconda ucsc-liftover
```

##  Running-the-pipeline
The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.

- Run pipeline directly by providing paths to the input files as arguments i.e
```
nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline -profile slurm -resume \
    --bloodCells '["baso", "rbc", "wbc"]' \
    --basePath "/new/path/to/data" \
    --ref19 "/new/path/to/hg37_snp_info_header.txt" \
    --ref38 "/new/path/to/hg38_snp_info_header.txt" \
    --chain_hg19_to_hg38 "/new/path/to/hg19ToHg38.over.chain.gz" \
    --chain_hg38_to_hg19 "/new/path/to/hg38ToHg19.over.chain.gz" \
    --target_genome_build 'hg19' \
    --plink_file '[["UGRC", "/new/path/to/uganda.bed", "/new/path/to/uganda.bim", "/new/path/to/uganda.fam"]]'
```

- Alternatively
Modify the conf/test.config file to suit the path to your data location, i.e
```
params.bloodCells = ["baso", "rbc", "wbc"]
params.basePath = "/new/path/to/data"
ref19 = "/new/path/to/hg37_snp_info_header.txt"
ref38 = "/new/path/to/hg38_snp_info_header.txt"
chain_hg19_to_hg38 = "/new/path/to/hg19ToHg38.over.chain.gz"
chain_hg38_to_hg19 = "/new/path/to/hg38ToHg19.over.chain.gz"
target_genome_build = 'hg38'
plink_file = [
    ['UGRC', '/new/path/to/uganda.bed', '/new/path/to/uganda.bim', '/new/path/to/uganda.fam']
]
```
 and then run the command as below:
   ```
   nextflow run nanjalaruth/PRS-Analysis-Nextflow-Pipeline -profile slurm -resume -c <path to your edited conf/test.config file> 
   ```

## To run the updated version of this pipeline, run:

 ```
 nextflow pull nanjalaruth/Intergrated_PRS_Analysis
 ```

## Workflow
A summary of the steps followed in our analysis include;
- Downloading score files from PGS catalogue
    - modify_score_file (Removes header from output file)
- Downloading metadata files from PGS catalogue
    - modify metadata file format   
    - modify_score_file_2 (Renaming rsID to chr_name:chr_position using the `snp info files` and metadata file)
- Liftover Genomic coordinates to match target genome build eg `hg38 to hg19`
    - modify_score_file_3 (modify format of liftover output)
- Computes PGS scores using `PLINK V1.9`
- Modify the output and concatenate scores for each phenotype
  
![pipeline](https://github.com/nanjalaruth/PRS-Analysis-Nextflow-Pipeline/blob/main/output/pipeline_info/output.png)

## Output
[Analysis_output](https://nanje.quarto.pub/intergrated_prs/)
[Nextflow report]()

## Support
I track open tasks using github's [issues](https://github.com/nanjalaruth/Intergrated_PRS_Analysis/issues)

## Citation
