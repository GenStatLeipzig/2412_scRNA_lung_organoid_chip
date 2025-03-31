#!/usr/bin/env Rscript
require(toolboxH)
require(here)


# prepare samples file
# Create a samplesheet CSV for nf-core/scrnaseq

# Get all fastq files recursively
allfiles <- dir("/net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip/",
    recursive = TRUE,
    pattern = "*.fastq.gz",
    full.names = TRUE
)

# Convert to data table
files_dt <- data.table(path = allfiles)

# Extract sample, read, and run information using stringr
files_dt[, basename := basename(path)]
files_dt[, sample := str_extract(basename, "^(S[1-4])")]
files_dt[, read := str_extract(basename, "R[12]")]
files_dt[, run := str_extract(basename, "_[001-2]+")]

files_dt[, .N, read]
files_dt[, .N, run]

# Create samplesheet according to nf-core/scrnaseq documentation
# Each sample will appear multiple times (for each run)
samplesheet <- dcast.data.table(
    files_dt,
    sample + run ~ read,
    value.var = "path"
)

# save
output_csv <- here("results/03_samplesheet.csv")

samplesheet_save <- samplesheet[, .(sample, fastq_1 = R1, fastq_2 = R2)]
samplesheet_save
# Write to CSV

fwrite(samplesheet_save, output_csv)
cat("Samplesheet created:", output_csv, "\n")
file.exists(output_csv)


# references
# https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads
# md5sum: a7b5b7ceefe10e435719edc1a8b8b2fa

# Check MD5 checksum of reference file
expected_md5 <- "a7b5b7ceefe10e435719edc1a8b8b2fa"
file_to_check <- here("results/cellranger1/refdata-gex-GRCh38-2024-A.tar.gz") # Replace with actual file path


actual_md5_pre <- system2("md5sum", file_to_check, stdout = TRUE)
actual_md5 <- substr(actual_md5_pre, 1, 32)
actual_md5 == expected_md5

# Clone the repository (if you haven't already)
# setwd(here("../nf_scrnaseq"))
# system("git clone https://github.com/nf-core/scrnaseq.git")



# teste pipeline

setwd(here("../nf_scrnaseq/scrnaseq"))

output_results <- here("results/cellranger1")
output_results

fasta_fn <- here("results/cellranger1/refdata-gex-GRCh38-2024-A/fasta/genome.fa")
file.exists(fasta_fn)

gtf_fn <- here("results/cellranger1/refdata-gex-GRCh38-2024-A/genes/genes.gtf")
file.exists(gtf_fn)



call2 <- paste(
    "nextflow run nf-core/scrnaseq --input ",
    output_csv,
    " --outdir ", output_results,
    " --aligner cellranger",
    " --protocol 10XV3",
    " --fasta", fasta_fn,
    " --gtf", gtf_fn,
    " -profile singularity"
)

system(call2)


finalizeSkript()

# knitr::opts_chunk$set(error = FALSE)
# filename_in = '/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/scripts/03_nf_scrnaseq.R'  # nolint
# file.exists(filename_in) # check if file exists
# knitr::spin(filename_in, knit = TRUE)
