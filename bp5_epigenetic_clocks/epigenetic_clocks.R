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

# One-time only step: You must download the PCClocks_data.qs2 file (~2GB).
# Creating the directory for storing the PCClocks file.
target_dir <- "path/to/save/methylCIPHER"
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# Downloading the PC file from Zenodo.  
options(timeout = 3600)   # allow enough time for the ~2 GB download
download.file(
	url = "https://zenodo.org/records/17162604/files/PCClocks_data.qs2?download=1",
	destfile = file.path(target_dir, "PCClocks_data.qs2"),
	mode = "wb"
)

# Read in the methylation file with Illumina cg id annotations
imputed_met_cpg <- read.table("/path/to/Illumina/ref/met.tsv", sep = "\t",
                                header = TRUE, row.names = 1)

# calcPCClocks() expects samples as rows, CpGs as columns
transposed <- t(imputed_met_cpg)

# Read in phenotype file. Columns must be formatted: Sample_ID, Female (1 = F, 0 = M), Age. 
pheno <- read.table("/path/to/pheno.tsv", sep = "\t", header = TRUE)

# Must load PC clock data before running calcPCClocks
pc_clocks_data <- load_PCClocks_data("/path/to/PCClocks_data.qs2")

# Estimate epigenetic age. Just put path to file holding PCC_clocks
result <- calcPCClocks(DNAm = transposed, pheno = pheno, RData = "/path/to/PCClocks_data.qs2")

# Save age estimates
write.table(result, "/path/to/epigenetic/age.tsv", sep = "\t",
quote = FALSE, row.names = TRUE, col.names = NA)

