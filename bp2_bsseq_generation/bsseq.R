# Installing libraries (Biostrings, BSseq, GenomicRanges, and BSgenome.Hsapiens.UCSC.hg38)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "Biostrings",
  "bsseq",
  "GenomicRanges",
  "BSgenome.Hsapiens.UCSC.hg38"
))

# Load Libraries 
library(Biostrings)
library(bsseq)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

# Restrict for autosomes
autosomes <- names(Hsapiens)[1:22]

# Find every CG dinucleotide position on each autosome.
cpg_positions_by_chr <- lapply(autosomes, function(chr) {
  start(matchPattern("CG", Hsapiens[[chr]]))
})

# Combine into one GRanges object covering every CpG site genome-wide (chr1-22).
all_cpg_sites <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cpg_positions_by_chr[[x]], width = 2))))

# Read in the Twist Targeted Methylome panel manifest (probe regions).
twist_manifest <- read.table(
  "/file/path",
  header = FALSE, sep = "\t",
  stringsAsFactors = FALSE
)[, 1:3]
colnames(twist_manifest) <- c("chr", "start", "end")

# Pad each probe by 200bp upstream/downstream. Captured DNA fragments are typically longer than the #120bp probe itself, so this recovers CpGs near probe edges that are still reliably captured.
twist_manifest_padded <- twist_manifest
twist_manifest_padded$start <- twist_manifest$start - 200
twist_manifest_padded$end   <- twist_manifest$end + 200

# Keep only genome-wide CpG sites that fall within the padded Twist regions. Sites outside these regions are off-target and unlikely to be high-quality.
cpgs_in_twist <- as.data.frame(
  GenomicRanges::intersect(x = all_cpg_sites, y = GRanges(twist_manifest_padded))
)
                                   
# width == 2 keeps only sites where the full CG dinucleotide overlaps the probe region (a width-1 overlap means only one base of the CG pair is captured, which we don't want).
cpgs_in_twist <- subset(cpgs_in_twist, width == 2)

# CX reports list the C and G of each CpG as separate lines (one per strand). Build one row for the plus-strand C position and one for the minus-strand G position, so both are present as loci when reading in the CX reports.
# (strandCollapse = TRUE, below, later merges these back into one CpG.)
plus_strand_sites <- cpgs_in_twist
plus_strand_sites$end <- plus_strand_sites$start
plus_strand_sites$strand <- "+"

minus_strand_sites <- plus_strand_sites
minus_strand_sites$start <- minus_strand_sites$end <- plus_strand_sites$end + 1
minus_strand_sites$strand <- "-"

cpg_loci <- rbind(plus_strand_sites, minus_strand_sites)
cpg_loci <- cpg_loci[order(cpg_loci$seqnames, cpg_loci$start), ]
cpg_loci$width <- 1

# Pull the CX reports together to generate a BSseq object
bismarkBSseq <- read.bismark(
  files = list.files("/path/to/CX/reports/", full.names = TRUE,
                      pattern = ".CX_report.txt.gz"),
  strandCollapse = TRUE,
  loci = GRanges(cpg_loci),
  verbose = TRUE,
  BACKEND = "HDF5Array",
  dir = "/path/to/Bsseq/object",
  replace = FALSE,   
  nThread = 8
)
