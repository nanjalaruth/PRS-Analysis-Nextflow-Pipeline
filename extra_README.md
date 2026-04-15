##  Data
2. SNP Info files for hg19/GRCh37 and hg38/GRCh38 genome builds
- The info files should be in a tsv format with header information, as shown in the examples below.
- If your file is in any other format, please convert it to this format using bash, python or R
- The examples are:
  - hg38 info files named as `hg38_snp_info_header.txt`

    | chr_name  | chr_position | rsID    |
    |-----------|--------------|---------|
    | 10   | 76684698    | rs241         |
    | 10   | 96480625    | rs243         |
    | 10   | 20703742    | rs244         |  

  - hg37 info files named as `hg37_snp_info_header.txt`
    | chr_name  | chr_position | rsID    |
    |-----------|--------------|---------|
    | 10   | 78444456    | rs241         |
    | 10   | 98240382     | rs243         |
    | 10   | 20992671    | rs244         |  

3. The liftover chain files
- Please download the chain files:
  - hg19ToHg38.over.chain.gz
  - hg38ToHg19.over.chain.gz

##  Tools

5. Liftover
We used conda for installation but feel free to install the liftover tool using other options
```
conda create -n liftover
conda activate liftover
conda install -c bioconda ucsc-liftover
```

##  Running-the-pipeline

### Required Arguments
| Argument  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| --basePath  | /new/path/to/all_blood_traits_prs_scores | Directory pattern for PGS Ids files      | 
| --bloodCells  | baso\rbc\wbc | Names of your phenotypes | 
| --ref19  | hg37_snp_info_header.txt | Genome build 37 SNP Info file with header information | 
| --ref38  | hg38_snp_info_header.txt | Genome build 38 SNP Info file with header information     | 
| --chain_hg19_to_hg38  | hg19ToHg38.over.chain.gz | hg19 to 38 chain file     | 
| --chain_hg38_to_hg19  | hg38ToHg19.over.chain.gz | hg38 to 19 chain file    | 
| --plink_file  | bed\bim\fam | PLINK genotype files   | 
| --target_genome_build| \<hg19\GRCh37>\,<hg38\GRCh38>| Path to the genome build to which the samples will be mapped |

- The pipeline does not require installation as `NextFlow` will automatically fetch it from `GitHub`.
- Modify the conf/test.config file particularly the lines below to suit the path to your data location:
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
