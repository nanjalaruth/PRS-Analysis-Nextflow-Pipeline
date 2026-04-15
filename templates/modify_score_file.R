#!/gpfs3/apps/eb/el8/2023a/skylake/software/R/4.3.2-gfbf-2023a/bin/Rscript

options(repos = c(CRAN = "https://cran.r-project.org"))

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr,
  data.table)

# Load the files
input_file <- read.table("${mod1}", header = TRUE, sep = "\\t")

#edit_file <- read.table("{edit_2}", header = TRUE, sep = "\\t")
#converted_file <- read.table("{edit_3}", header = TRUE, sep = "\\t")

# Check if the files are similar
#are_files_similar <- identical(edit_file, converted_file)

#if (are_files_similar) {
  # If the files are identical, just use one of them (e.g., edit_2)
  #print("Files are identical, using edit_2 directly.")
  #final_result <- edit_file
#} else {
  # If the files are not identical, proceed with the merge
  #print("Files are not identical, proceeding with merge.")
  
  # Merge the files based on the SNP column
  #merged_file <- merge(edit_file, converted_file, by = "SNP", all.x = TRUE)
  
  # Rename SNP column to Old_SNP
  #merged_file <- merged_file %>%
    #dplyr::rename(Old_SNP = SNP)
  
  # Create a new SNP column by combining lifted_chr_name and end_position without the "chr" prefix
  #merged_file\$SNP <- paste(merged_file\$chr_name, merged_file\$end_position, sep = ":")
  
  # Use merged result as final output
  #final_result <- merged_file
#}

# create a snp column
if (!"SNP" %in% names(input_file)) {
  chr_name_index <- which(names(input_file) == "hm_chr")
  chr_pos_index  <- which(names(input_file) == "hm_pos")
  
  if (length(chr_name_index) == 0 || length(chr_pos_index) == 0) {
    stop("Columns 'chr_name' or 'chr_position' not found.")
  }
  
  input_file\$SNP <- paste0(input_file[[chr_name_index]], ":", input_file[[chr_pos_index]])
}

scoring_file <- data.frame(SNP = input_file\$SNP, A1 = input_file\$effect_allele, BETA = input_file\$effect_weight)
scoring_file_uniq <- scoring_file[!duplicated(scoring_file\$SNP), ]
score_file <- na.omit(scoring_file_uniq)
clean_score_file <- score_file[score_file[, 1] != "", , drop = FALSE]
write.table(clean_score_file, file = "${blood_trait}_${pgs_id}_modified.txt", sep = "\\t", col.names = T, row.names = F, quote = F)

