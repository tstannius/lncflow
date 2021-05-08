#!/usr/bin/env Rscript

# usage
#   this script converts a loom file to an rds file

# based on the Seurat conversion vignette: 
#   https://satijalab.org/seurat/articles/conversion_vignette.html

# PARSE ARGUMENTS #############################################################
suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option(c("-i", "--input"), type="character",
                  help="Input .loom file."),
    make_option(c("-o", "--output"), type="character",
                  help="Output .rds file.")
)

args <- parse_args(OptionParser(option_list=option.list))

if (is.null(args$input)) {
  stop("[ERROR]: No input .loom file provided. See `--help`.")
}
if (is.null(args$output)) {
  stop("[ERROR]: No output .rds file path provided. See `--help`.")
}

# e.g. "/scratch/gaurav/velocyto_runs/d14v/d14v-5000_cells/velocyto/d14v-5000_cells.loom"
path.loom <- args$input
path.rds <- args$output
project_id <- tools::file_path_sans_ext(basename(path.loom))

sprintf("in: %s", path.loom)
sprintf("out: %s", path.rds)

# RUN WORKFLOW ################################################################
# if optparse completes, continue with rest of script
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk))
set.seed(42)


data.loom <- SeuratDisk::Connect(filename = path.loom, mode = "r")
data.rds <- as.Seurat(data.loom)
saveRDS(data.rds, file = path.rds)
data.loom$close_all()

sessionInfo()
# savehistory(paste0("out/", project_id, "_sessioninfo.txt"))
