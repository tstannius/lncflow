#!/usr/bin/env Rscript

# usage
#   this script runs Seurat::FindAllMarkers() on the input Seurat Object rds

# based on the following vignettes and guides:
#   guided clustering:      https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# notes
# - none atm

# PARSE ARGUMENTS #############################################################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))

option.list <- list(
    make_option(c("-i", "--input_rds"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("-o", "--output"), type="character",
                  help="Output .csv filename.")
)

args <- parse_args(OptionParser(option_list=option.list))

if (is.null(args$input_rds)) {
  stop("[ERROR]: No input .rds file provided. See `--help`.")
}
if (is.null(args$output)) {
  stop("[ERROR]: No output filename provided. See `--help`.")
}


path.rds.in <- args$input
path.csv.out <- args$output
# project.id <- tools::file_path_sans_ext(basename(path.rds.in))

# RUN WORKFLOW ################################################################
source("/nfsdata/projects/tstannius/lncflow/scripts/seurat/helpers.R")
print.status("[INFO]: Loading dependencies.")
suppressPackageStartupMessages(library(Seurat))

print.status(paste("[INFO]: Reading in:", path.rds.in))
sobj <- readRDS(path.rds.in)
DefaultAssay(sobj) <- "RNA"

print.status("[INFO]: Running FindAllMarkers()")
markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print.status(paste("[INFO]: Saving markers to:", path.rds.in))
write.csv(markers, path.csv.out)
