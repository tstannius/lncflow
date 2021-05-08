#!/usr/bin/env Rscript

# usage
# this script converts a loom file to an rds file

# based on the Seurat conversion vignette: https://satijalab.org/seurat/v3.1/conversion_vignette.html

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
  stop("Input file not specified. Usage: loom2rds.R path/to/data.loom", call.=FALSE)
}

path_loom <- args[1]
project_id <- tools::file_path_sans_ext(basename(path_loom))
path_rds <- paste0("out/", project_id, "_raw.rds")

dir.create(file.path(getwd(), "out"), showWarnings = FALSE)

sprintf("in: %s", path_loom)
sprintf("out: %s", path_rds)

# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)


data.loom <- connect(filename = path_loom, mode = "r+", skip.validate = TRUE)
data.rds <- as.Seurat(data.loom)
saveRDS(data.rds, file = path_rds)
data.loom$close_all()

sessionInfo()
savehistory(paste0("out/", project_id, "_sessioninfo.txt"))
