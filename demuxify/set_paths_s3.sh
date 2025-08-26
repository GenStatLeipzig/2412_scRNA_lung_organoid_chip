#!/bin/bash
# filepath: demuxify/set_paths_s3.sh
# source set_paths_s3.sh before running demuxify/Demuxafy
# This script sets the paths for the input files and output directories for the demuxify/Demuxafy pipeline
# Test if variables are set
# echo "JSON: $JSON"
# echo "N_DOUB: $N_DOUB"
# source demuxify/set_paths_s3.sh && mkdir -p "$OUTDIR" "$SCDS_OUTDIR" "$SCDBLFINDER_OUTDIR" "$SCRUBLET_OUTDIR" "$DOUBLETFINDER_OUTDIR" "$DOUBLETDETECTION_OUTDIR" "$SOLO_OUTDIR" 

# source demuxify/set_paths_s3.sh && singularity exec --bind /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/ demuxify/Demuxafy.sif solo -o $SOLO_OUTDIR -j $JSON -d $COUNTS -e $N_DOUB  2>&1 | tee $SOLO_OUTDIR/run.log

# singularity exec --bind /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/ demuxify/Demuxafy.sif solo_summary.py -b $COUNTS/barcodes.tsv.gz -s $SOLO_OUTDIR


# source demuxify/set_paths_s3.sh && singularity exec --bind /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/ demuxify/Demuxafy.sif Scrublet.py -m $COUNTS -o $SCRUBLET_OUTDIR -f $BARCODES 2>&1 | tee $SCRUBLET_OUTDIR/run.log

# source demuxify/set_paths_s3.sh && singularity exec --bind /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/ demuxify/Demuxafy.sif scDblFinder.R -o $SCDBLFINDER_OUTDIR -t $COUNTS -b $BARCODES_FILTERED 2>&1 | tee $SCDBLFINDER_OUTDIR/run.log

# source demuxify/set_paths_s3.sh && singularity exec --bind /net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/ demuxify/Demuxafy.sif scds.R -o $SCDS_OUTDIR -t $COUNTS -b $BARCODES 2>&1 | tee $SCDS_OUTDIR/run.log


# VCF=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/TestData4PipelineFull/test_dataset.vcf # only needed for genetics-based methods

COUNTS=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/results/cellranger1/cellranger/count/S3/outs/filtered_feature_bc_matrix ## Change this based on the path on your system
BARCODES=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/results/cellranger1/cellranger/S3/cellbender_removebackground/S3_cell_barcodes.csv ## Change this based on the path on your system
BARCODES_FILTERED=$BARCODES 

# BAM=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/TestData4PipelineFull/test_dataset/possorted_genome_bam.bam # only needed for genetics-based methods
# INDS=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/TestData4PipelineFull/donor_list.txt # only needed for genetics-based methods
# GTF=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/TestData4PipelineFull/genes.gtf # only needed for genetics-based methods
DOUBLETS=800 # Change this based on the number of doublets you expect in your dataset,used from multiQC report from cellranger entered in https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/test.html
N_DOUB=$DOUBLETS
OUTDIR=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/S3
# DEMUXALOT_OUTDIR=$OUTDIR/demuxalot # only needed for genetics-based methods
# DROPULATION_OUTDIR=$OUTDIR/dropulation # only needed for genetics-based methods
SCDS_OUTDIR=$OUTDIR/scds
SCDBLFINDER_OUTDIR=$OUTDIR/scDblFinder
SCRUBLET_OUTDIR=$OUTDIR/scrublet
DOUBLETFINDER_OUTDIR=$OUTDIR/doubletFinder
DOUBLETDETECTION_OUTDIR=$OUTDIR/DoubletDetection
SOLO_OUTDIR=$OUTDIR/solo
JSON=/net/ifs1/san_projekte/projekte/genstat/02_projekte/2412_scRNA_lung_organoid_chip/gitOrganoid_hk/demuxify/solo_model.json