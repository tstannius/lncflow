#!/usr/bin/env Rscript

# usage
#   this script demultiplexes a Seurat object

# based on the Seurat Demultiplexing vignette:
#   https://satijalab.org/seurat/articles/hashing_vignette.html

# PARSE ARGUMENTS #############################################################
suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option(c("-r", "--seurat_rds"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("-c", "--cite_seq_dir"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("-o", "--output"), type="character",
                  help="Output .rds file.")
)

args <- parse_args(OptionParser(option_list=option.list))

if (is.null(args$seurat_rds)) {
  stop("[ERROR]: No input .rds file provided. See `--help`.")
}
if (is.null(args$cite_seq_dir)) {
  stop("[ERROR]: No input CITE-seq directory provided. See `--help`.")
}
if (is.null(args$output)) {
  stop("[ERROR]: No output .rds file path provided. See `--help`.")
}


path.rds.in <- args$seurat_rds
path.citeseq.dir <- args$cite_seq_dir
path.rds.out <- args$output
project_id <- tools::file_path_sans_ext(basename(path.rds.in))

# RUN WORKFLOW ################################################################
# if optparse completes, continue with rest of script
suppressPackageStartupMessages(library(Seurat))
set.seed(42)

# Load in the UMI matrix
umis <- readRDS(path.rds.in)

# Read in HTO data
# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
htos <- Read10X(path.citeseq.dir, gene.column = 1)

# process HTO data
htos <- htos[!(rownames(htos) %in% c("unmapped")),] # Throw away any unmapped rows

prefix.dset <- unlist(strsplit(colnames(umis)[1], ':'))[1]
colnames(htos) <- paste(prefix.dset, ':', colnames(htos),'x', sep='') # to match umis

# Determine joint names between umis and hashtag
joint.bcs <- intersect(colnames(umis), colnames(htos))

# Subset RNA and HTO counts by joint cell barcodes
umis <- subset(umis, cells = joint.bcs)
htos <- as.matrix(htos[, joint.bcs])

# Setup Seurat object for running Seurat::HTODemux()
hashtag <- umis

# Add HTO data as a new assay independent from RNA
hashtag[['HTO']] <- CreateAssayObject(counts = htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
hashtag <- NormalizeData(hashtag, assay = 'HTO', normalization.method = 'CLR')

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag, assay = 'HTO', positive.quantile = 0.99)

# Global classification results
table(hashtag$HTO_classification.global)

# Extract singlets, i.e. those values that unambiguously belong to one cell
Idents(hashtag) <- "HTO_classification.global" # ensure that correct idents is selected
singlet <- subset(hashtag, idents = "Singlet")

saveRDS(singlet, path.rds.out)
