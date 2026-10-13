# Install EpiDISH if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
if (!requireNamespace("EpiDISH", quietly = TRUE)) {
BiocManager::install("EpiDISH")
}

library(EpiDISH)
library(data.table)

# Load the EpiDISH blood reference matrix. Here we use centDHSbloodDMC.
data(centDHSbloodDMC.m)

# Read in the methylation matrix (CpG ID row names if using the Illumina/EPIC reference
met_df <- read.table("/path/to/Illumina/ref/met.tsv", sep = "\t",
header = TRUE, row.names = 1)
met_matrix <- as.matrix(met_df)

# Estimate cell type proportions via robust partial correlation (RPC).
BloodFrac.m_output <- epidish(beta.m = met_matrix,
ref.m = centDHSbloodDMC.m,	# or other methylation reference
method = "RPC")$estF

# Saving cell proportion 
cell_props_df <- as.data.frame(BloodFrac.m_output)
fwrite(cell_props_df, "/path/to/cell/fractions.tsv", sep = "\t", row.names = TRUE)

