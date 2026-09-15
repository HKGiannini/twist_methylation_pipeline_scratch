# Install methylCIPHER package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("methylCIPHER", quietly = TRUE)) {
  devtools::install_github("HigginsChenLab/methylCIPHER")
}

# Loading libraries
library(methylCIPHER)
library(tibble)
library(dplyr)

# One-time only step: You must download the PCClocks_data.qs2 file (~2GB). This file contains training data for the clocks. 
# Creating the directory for storing the PCClocks file.
target_dir <- get_methylCIPHER_path()
target_dir  # print it so you can see the path
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# Downloading the file from Zenodo.
download.file(
  url = "https://zenodo.org/records/17162604/files/PCClocks_data.qs2?download=1",
  destfile = file.path(target_dir, "PCClocks_data.qs2"),
  mode = "wb"
)

# Read in the beta file Illumina Cpg annotated (see BP3) since PCClocks was built using the EPIC array. 
imputed_beta_cpg <- read.table("/path/to/Illumina/ref/beta.tsv", sep = "\t",
                                header = TRUE, row.names = 1)

# calcPCClocks() expects samples as rows, CpGs as columns. You will need to transpose the beta file.
transposed <- t(imputed_beta_cpg)

# Read in phenotype file. It must be in the format: Sample_ID, Female (1 = Female, 0 = Male), Age. You need this info because some PC clocks (e.g. PCGrimAge) use chronological age and/or sex as model inputs, not just methylation values.
pheno <- read.table("/path/to/pheno.tsv", sep = "\t", header = TRUE)

# Estimate epigenetic age 
result <- calcPCClocks(DNAm = transposed, pheno = pheno, RData = "/path/to/PCClocks_data.qs2")

# Save age estimates
write.table(result, "/path/to/epigenetic/age.tsv", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)





