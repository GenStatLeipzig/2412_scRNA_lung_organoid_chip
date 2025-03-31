

``` r
library(toolboxH)
```

```
## Lade nötiges Paket: data.table
```

```
## data.table 1.16.0 using 64 threads (see ?getDTthreads).  Latest news: r-datatable.com
## Lade nötiges Paket: R.utils
## Lade nötiges Paket: R.oo
## Lade nötiges Paket: R.methodsS3
## R.methodsS3 v1.8.2 (2022-06-13 22:00:14 UTC) successfully loaded. See ?R.methodsS3 for help.
## R.oo v1.26.0 (2024-01-24 05:12:50 UTC) successfully loaded. See ?R.oo for help.
## 
## Attache Paket: 'R.oo'
## 
## Das folgende Objekt ist maskiert 'package:R.methodsS3':
## 
##     throw
## 
## Die folgenden Objekte sind maskiert von 'package:methods':
## 
##     getClasses, getMethods
## 
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     attach, detach, load, save
## 
## R.utils v2.12.3 (2023-11-18 01:00:02 UTC) successfully loaded. See ?R.utils for help.
## 
## Attache Paket: 'R.utils'
## 
## Das folgende Objekt ist maskiert 'package:utils':
## 
##     timestamp
## 
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     cat, commandArgs, getOption, isOpen, nullfile, parse, use, warnings
## 
## Lade nötiges Paket: fdrtool
## Lade nötiges Paket: png
## Lade nötiges Paket: RColorBrewer
## Lade nötiges Paket: readxl
## Lade nötiges Paket: scales
## Lade nötiges Paket: stringr
## Lade nötiges Paket: testthat
## 
## Attache Paket: 'testthat'
## 
## Das folgende Objekt ist maskiert 'package:R.oo':
## 
##     equals
## 
## Lade nötiges Paket: eulerr
## 
## Attache Paket: 'toolboxH'
## 
## Das folgende Objekt ist maskiert 'package:eulerr':
## 
##     venn
```

``` r
library(here)
```

```
## here() starts at /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk
```

``` r
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
```

```
## Concatenating files forS1R1->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S1_S23_R1_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S1/S1_S23_R1_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S1/S1_S23_R1_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S1_S23_R1_001.fastq.gz
```

```
## Concatenating files forS1R2->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S1_S23_R2_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S1/S1_S23_R2_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S1/S1_S23_R2_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S1_S23_R2_001.fastq.gz
```

```
## Concatenating files forS2R1->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S2_S24_R1_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S2/S2_S24_R1_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S2/S2_S24_R1_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S2_S24_R1_001.fastq.gz
```

```
## Concatenating files forS2R2->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S2_S24_R2_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S2/S2_S24_R2_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S2/S2_S24_R2_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S2_S24_R2_001.fastq.gz
```

```
## Concatenating files forS3R1->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S3_S25_R1_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S3/S3_S25_R1_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S3/S3_S25_R1_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S3_S25_R1_001.fastq.gz
```

```
## Concatenating files forS3R2->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S3_S25_R2_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S3/S3_S25_R2_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S3/S3_S25_R2_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S3_S25_R2_001.fastq.gz
```

```
## Concatenating files forS4R1->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S4_S26_R1_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S4/S4_S26_R1_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S4/S4_S26_R1_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S4_S26_R1_001.fastq.gz
```

```
## Concatenating files forS4R2->/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S4_S26_R2_001.fastq.gz
## Using command:
```

```
## cat /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S4/S4_S26_R2_001.fastq.gz /net/ifs1/san_projekte/projekte/genstat/01_daten/2412_scRNA_lung_organoid_chip//Karen_fastq_bothruns/S4/S4_S26_R2_002.fastq.gz > /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined/S4_S26_R2_001.fastq.gz
```

``` r
message("All files have been concatenated and saved to", output_dir, "\n")
```

```
## All files have been concatenated and saved to/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/../fastq_combined
```

``` r
finalizeSkript()
```

```
## ==================================================================================
## 
## 
## Warnings found so far:
```

```
## 
## kann Datei '/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/scripts/02_concatenate_fastq.R' nicht öffnen: Datei oder Verzeichnis nicht gefunden 
##                                                                                                                                                                                                1
```

```
## ==================================================================================
## 
## 
## Session Info::
```

```
## R version 4.4.2 (2024-10-31)
## Platform: x86_64-suse-linux-gnu
## Running under: openSUSE Leap 15.6
## 
## Matrix products: default
## BLAS:   /usr/lib64/R/lib/libRblas.so 
## LAPACK: /usr/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Berlin
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] here_1.0.1         toolboxH_0.2.17    eulerr_7.0.2       testthat_3.2.1.1  
##  [5] stringr_1.5.1      scales_1.3.0       readxl_1.4.3       RColorBrewer_1.1-3
##  [9] png_0.1-8          fdrtool_1.2.18     R.utils_2.12.3     R.oo_1.26.0       
## [13] R.methodsS3_1.8.2  data.table_1.16.0 
## 
## loaded via a namespace (and not attached):
##  [1] vctrs_0.6.5      cli_3.6.3        knitr_1.48       rlang_1.1.4     
##  [5] xfun_0.46        stringi_1.8.4    glue_1.7.0       colorspace_2.1-1
##  [9] rprojroot_2.0.4  brio_1.1.5       grid_4.4.2       cellranger_1.1.0
## [13] evaluate_0.24.0  munsell_0.5.1    lifecycle_1.0.4  compiler_4.4.2  
## [17] Rcpp_1.0.13      R6_2.5.1         magrittr_2.0.3   tools_4.4.2
```

```
## ==================================================================================
## 
## 
## Total Time:
```

``` r
# knitr::opts_chunk$set(error = FALSE)
# filename_in = '/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/scripts/02_concatenate_fastq.R'  # nolint
# file.exists(filename_in) # check if file exists
# knitr::spin(filename_in, knit = TRUE)
```

