library(toolboxH)
library(here)


# get list of files
allfiles <- dir("/net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip/", recursive = T, pattern = "*.fastq.gz", full.names = TRUE)



# Convert the file paths to a data.table
files_dt <- data.table(path = allfiles)

# Extract sample, read, and run information
files_dt[, c("sample", "read", "run") := {
    # Extract basename and split components
    base <- basename(path)
    # Extract sample ID (S1, S2, etc.)
    sample <- gsub("(S[1-4])_.*", "\\1", base)
    # Extract read type (R1 or R2)
    read <- gsub(".*_(R[12])_.*", "\\1", base)
    # Extract run number (001 or 002)
    run <- gsub(".*_(00[12]).fastq.gz", "\\1", base)
    list(sample, read, run)
}]

# Create output directory
output_dir <- here("../fastq_combined")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Group by sample and read to identify pairs to concatenate
sample_read_groups <- files_dt[, .(paths = list(path)), by = .(sample, read)]

# Concatenate each pair
for (i in 1:nrow(sample_read_groups)) {
    # Get files for this sample and read
    #   i=1
    current_files <- sample_read_groups[i, paths][[1]]

    # Sort files to ensure 001 comes before 002
    current_files <- sort(current_files)

    # Create output filename (using 001 convention)
    output_base <- basename(current_files[1])
    output_file <- file.path(output_dir, output_base)

    # Concatenate files
    message(
        "Concatenating files for", sample_read_groups[i, sample],
        sample_read_groups[i, read], "->", output_file, "\nUsing command:\n"
    )

    # Use cat command to concatenate
    cmd <- paste("cat", paste(current_files, collapse = " "), ">", output_file)
    cat(cmd)
    system(cmd)
}

message("All files have been concatenated and saved to", output_dir, "\n")


finalizeSkript()

# knitr::opts_chunk$set(error = FALSE)
# filename_in = '/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/scripts/02_concatenate_fastq.R'  # nolint
# file.exists(filename_in) # check if file exists
# knitr::spin(filename_in, knit = TRUE)
