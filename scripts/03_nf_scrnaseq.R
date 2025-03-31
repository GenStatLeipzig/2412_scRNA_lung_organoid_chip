require(toolboxH)
require(here)
# prepare samples file
# Create a samplesheet CSV for nf-core/scrnaseq
library(data.table)

combined_dir <- here("../fastq_combined")
output_csv <- paste0(combined_dir, "/samplesheet.csv")

# Get list of FASTQ files
fastq_files <- list.files(path = combined_dir, pattern = "*.fastq.gz", full.names = TRUE)

# Extract sample information
samples <- data.table(
    fastq_path = fastq_files,
    sample = str_extract(basename(fastq_files), "^[^_]+")
)

# Create unique sample entries
unique_samples <- unique(samples$sample)
unique_samples

# Create samplesheet

samplesheet <- lapply(unique_samples, function(s) {
    # s = unique_samples[1]
    sample_files <- samples[sample == s, fastq_path]

    # Identify read types (R1, R2, I1)
    r1_file <- sample_files[grep("_R1_", sample_files)]
    r2_file <- sample_files[grep("_R2_", sample_files)]
    # i1_file <- sample_files[grep("_I1_", sample_files)] # already demultiplexed, i.e. no multiple libraries from different project from seq center present

    data.table(
        sample = s,
        fastq_1 = r1_file,
        fastq_2 = r2_file
        # ,
        # expected_cells = 10000 # Note that since cellranger v7, it is not recommended anymore to supply the --expected-cells parameter.
        # ,
        # fastq_3 = i1_file
    )
}) %>% rbindlist()
samplesheet
# Write to CSV
fwrite(samplesheet, output_csv)
cat("Sample sheet created:", output_csv, "\n")

# Clone the repository (if you haven't already)
setwd(here("../nf_scrnaseq"))
system("git clone https://github.com/nf-core/scrnaseq.git")

# test the pipeline

setwd(here("../nf_scrnaseq/scrnaseq"))
# system("nextflow run nf-core/scrnaseq -profile test,singularity --outdir minitest --validate_params=false")
# does not run, error
# RuntimeError: cannot cache function 'sparse_mean_var_minor_axis': no locator available for file '/usr/local/lib/python3.8/site-packages/scanpy/preprocessing/_utils.py'


# Run the pipeline with your combined FASTQ files
# nextflow run main.nf \
#   --input samplesheet.csv \
#   --genome GRCh38 \
#   --phenotype_data <path_to_phenotype_file> \
#   --outdir <output_directory> \
#   --protocol 10XV3 \
#   -profile docker

finalizeSkript()

# knitr::opts_chunk$set(error = FALSE)
# filename_in = '/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/scripts/03_nf_scrnaseq.R'  # nolint
# file.exists(filename_in) # check if file exists
# knitr::spin(filename_in, knit = TRUE)
