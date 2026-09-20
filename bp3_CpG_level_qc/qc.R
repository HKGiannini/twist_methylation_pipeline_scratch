# Install comethyl (GitHub-hosted, not on Bioconductor/CRAN) if not already installed
if (!requireNamespace("comethyl", quietly = TRUE)) {
  BiocManager::install("cemordaunt/comethyl")
}

library(comethyl)   # filterCpGs()
library(dplyr)      
library(bsseq)
library(HDF5Array)  

# Load the HDF5-backed BSseq object generated in BP2
bsseq_obj <- loadHDF5SummarizedExperiment(dir = "/path/to/Bsseq/object")

# Keep only CpGs with >5x coverage in at least 75% of samples (PMID: 40402999)
# Adjust `cov`/`perSample` per your downstream analysis plan (see BP3 intro).
filtered <- filterCpGs(
  bsseq_obj,
  cov = 5,
  perSample = 0.75,
  file = "/path/to/filtered/Bsseq/object"
)

# Extract genomic coordinates from the filtered object
positions <- granges(filtered)
chrom <- as.character(seqnames(positions))
pos   <- start(positions)

# Extract raw methylation levels and coverage
methylation_df <- data.frame(chromosome = chrom, position = pos,
                              as.data.frame(getMeth(filtered, type = "raw")))
coverage_df    <- data.frame(chromosome = chrom, position = pos,
                              as.data.frame(getCoverage(filtered, type = "Cov")))

sample_cols <- setdiff(colnames(methylation_df), c("chromosome", "position"))

# Impute missing methylation values to the per-CpG (row) mean across samples.
methylation_df_imputed <- methylation_df
cpg_means <- rowMeans(methylation_df[sample_cols], na.rm = TRUE)

for (col in sample_cols) {
  na_idx <- is.na(methylation_df_imputed[[col]])
  methylation_df_imputed[[col]][na_idx] <- cpg_means[na_idx]
}

# Confirm no NAs remain
sum(is.na(methylation_df_imputed[sample_cols]))

# Read in Illumina's manifest. 
zip_path <- "/path/to/Illumina/manifest/zip"
csv_name <- unzip(zip_path, list = TRUE)$Name[1]
unzip(zip_path, exdir = "data/manifest")

# skip = 7: the manifest file has 7 lines of header metadata before the real column headers
manifest <- read.csv(file.path("data/manifest", csv_name),
                      skip = 7, stringsAsFactors = FALSE)

# Standardize chromosome naming (ensure "chr" prefix) and apply the +1 offset to align Illumina's coordinate convention with hg38.
manifest <- manifest %>%
  mutate(CHR_hg38 = ifelse(startsWith(CHR_hg38, "chr"), CHR_hg38, paste0("chr", CHR_hg38)),
         Start_hg38 = as.numeric(Start_hg38) + 1)

# Map chr:pos to Illumina CpG IDs (cgXXXXXXXX) for both methylation and coverage matrix
filtered_met_imputed_cpgid <- methylation_df_imputed %>%
  inner_join(manifest %>% select(IlmnID, CHR_hg38, Start_hg38),
             by = c("chromosome" = "CHR_hg38", "position" = "Start_hg38")) %>%
  relocate(IlmnID)

filtered_coverage_cpgid <- coverage_df %>%
  inner_join(manifest %>% select(IlmnID, CHR_hg38, Start_hg38),
             by = c("chromosome" = "CHR_hg38", "position" = "Start_hg38")) %>%
  relocate(IlmnID)

# Save the met and coverage matrices as both the Twist (chr:pos) and Illumina (cpgXXXXXXXX) formats. Will give you 4 different files. 
coverage_df_out <- coverage_df %>%
  mutate(`chr:pos` = paste0(chromosome, ":", position)) %>%
  select(`chr:pos`, everything(), -chromosome, -position)

methylation_df_imputed_out <- methylation_df_imputed %>%
  mutate(`chr:pos` = paste0(chromosome, ":", position)) %>%
  select(`chr:pos`, everything(), -chromosome, -position)

filtered_coverage_cpgid_out <- filtered_coverage_cpgid %>%
  select(-chromosome, -position)

filtered_met_imputed_cpgid_out <- filtered_met_imputed_cpgid %>%
  select(-chromosome, -position)

# Coverage, Twist chr:pos
write.table(coverage_df_out, "path/to/twist/ref/coverage.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(coverage_df_out, "path/to/twist/ref/coverage.rds")

# Coverage, CpG ID only
write.table(filtered_coverage_cpgid_out, "path/to/Illumina/ref/coverage.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(filtered_coverage_cpgid_out, "path/to/Illumina/ref/coverage.rds")

# Imputed mets, Twist chr:pos
write.table(methylation_df_imputed_out, "path/to/twist/ref/met.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(methylation_df_imputed_out, "path/to/twist/ref/met.rds")

# Imputed mets, CpG ID only
write.table(filtered_met_imputed_cpgid_out, "path/to/Illumina/ref/met.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(filtered_met_imputed_cpgid_out, "path/to/Illumina/ref/met.rds")
