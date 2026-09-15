# Core R dependencies

R version: **4.4.2** (2024-10-31), Ubuntu 24.04.1 LTS

This is a trimmed list of the packages the pipeline directly depends on (BP2–BP6). The full dependency tree, including everything pulled in transitively, is recorded verbatim in `sessionInfo.txt` for exact reproducibility.

| Package | Version | Used in |
|---|---|---|
| bsseq | 1.42.0 | BP2 |
| GenomicRanges | 1.58.0 | BP2 |
| Biostrings | 2.74.1 | BP2 |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 | BP2 |
| BSgenome | 1.74.0 | BP2 (base package) |
| data.table | 1.18.4 | BP3 (matrix formatting) |
| Matrix | 1.7-1 | BP3 |
| HDF5Array / rhdf5 | 1.34.0 / 2.50.2 | BP3 (large matrix storage) |
| SummarizedExperiment | 1.36.0 | BP3 |
| EpiDISH | 2.22.0 | BP4 |
| methylCIPHER | 0.2.0 | BP5 |
| dplyr | 1.2.1 | general |
