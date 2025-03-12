#!/gpfs3/apps/eb/el8/2023a/skylake/software/R/4.3.2-gfbf-2023a/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr,
  data.table)

# Get list of file names ending with .tsv in the directory
tsv_files <- list.files(pattern = "\\\\_pgs_scores.txt\$")

# Read and merge the TSV files
merged_data <- lapply(tsv_files, function(file) read.table(file, header = TRUE, sep = "\\t")) %>%
  Reduce(function(x, y) merge(x, y, by = "FID"), .)

# Write combined data to a new TSV file
write.table(merged_data, file = "${dataset}_${blood_trait}_pgs_scores.txt", sep = "\\t", row.names = FALSE, quote = FALSE)
