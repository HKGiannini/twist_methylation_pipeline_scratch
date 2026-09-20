# Install EpiDISH if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("EpiDISH", quietly = TRUE)) {
  BiocManager::install("EpiDISH")
}

# Load libraries
library(EpiDISH)
library(data.table)

# Read in the methylation matrix (CpG ID row names if using the Illumina/EPIC reference; chr:pos if using the Twist reference — see BP3).
met_df <- read.table("/path/to/Illumina/ref/met.tsv", sep = "\t",
                       header = TRUE, row.names = 1)
met_matrix <- as.matrix(met_df)

# Estimate cell type proportions via robust partial correlation (RPC).
BloodFrac.m_output <- epidish(met.m = met_matrix,
                               ref.m = X,
                               method = "RPC")$estF

# Save cell props 
cell_props_df <- as.data.frame(BloodFrac.m_output)
fwrite(cell_props_df, "/path/to/cell/props.tsv", sep = "\t", row.names = TRUE)
