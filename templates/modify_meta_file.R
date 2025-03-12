#!/gpfs3/apps/eb/el8/2023a/skylake/software/R/4.3.2-gfbf-2023a/bin/Rscript

options(repos = c(CRAN = "https://cran.r-project.org"))

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr,
  data.table)

result = read.csv("${meta_files}", header = TRUE)

scoring_file <- data.frame( blood_trait = result\$Reported.Trait,  pgs_id = result\$Polygenic.Score..PGS..ID, Genome_Build = result\$Original.Genome.Build)
#scoring_file_uniq <- scoring_file[!duplicated(scoring_file\$SNP), ]
#score_file <- na.omit(scoring_file_uniq)
#clean_score_file <- score_file[score_file[, 1] != "", , drop = FALSE]
write.table(scoring_file, file = "${blood_trait}_${pgs_id}_meta_modified.txt", sep = "\\t", col.names = T, row.names = F, quote = F)

