#!/usr/bin/env Rscript

# usage
#   this script runs Batchelor FastMNN to integrate two datasets.
#   it then identifiers marker genes for each dataset
#   finally it saves the markers for either dataset and the integrated dataset

# PARSE ARGUMENTS #############################################################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))

option.list <- list(
    make_option(c("--input_rds_1"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("--input_rds_2"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("--output_csv_1"), type="character",
                  help="Output .csv filename."),
    make_option(c("--output_csv_2"), type="character",
                  help="Output .csv filename."),
    make_option(c("--output_rds_integrated"), type="character",
                  help="Output .rds integrated filename.")
)

args <- parse_args(OptionParser(option_list=option.list))

if (is.null(args$input_rds_1) | is.null(args$input_rds_2)) {
  stop("[ERROR]: Please specify 2 input .rds files. See `--help`.")
}
if (is.null(args$output_csv_1) | is.null(args$output_csv_2) | is.null(args$output_rds_integrated)) {
  stop("[ERROR]: Please specify output_csv and output_rds. See `--help`.")
}


# SETUP VARIABLES #############################################################
path.rds.1 <- args$input_rds_1 # e.g. "/nfsdata/projects/tstannius/lncflow/out/hg38/d14D/d14D-nomux-seurat_preproc.rds"
path.rds.2 <- args$input_rds_2 # e.g. "/nfsdata/projects/tstannius/lncflow/out/hg38/d14V/d14V-demux-seurat_preproc.rds"

out.csv.1 <- args$output_csv_1 # e.g. "/nfsdata/projects/tstannius/lncflow/out/hg38/d14D/d14D-nomux-markers_region.csv"
out.csv.2 <- args$output_csv_2 # e.g. "/nfsdata/projects/tstannius/lncflow/out/hg38/d14V/d14V-demux-markers_region.csv"

out.rds <- args$output_rds_integrated # e.g. "/nfsdata/projects/tstannius/lncflow/out/hg38/d14/d14-integrated_fastMNN.rds"

project.id.1 <- unlist(str_split(tools::file_path_sans_ext(basename(path.rds.1)), '-'))
region.id.1 <- project.id.1[1]
project.id.2 <- unlist(str_split(tools::file_path_sans_ext(basename(path.rds.2)), '-'))
region.id.2 <- project.id.2[1]
integrated.id <- basename(dirname(out.rds))

p.list <- list() # list for plots


# RUN WORKFLOW ################################################################
source("/nfsdata/projects/tstannius/lncflow/scripts/seurat/helpers.R")
print.status("[INFO]: Loading dependencies.")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SeuratWrappers)) # contains the Seurat Batchelor wrapper for fastMNN
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# Read data ###################################################################
# N.B. the day 14 dorsal dataset was previously integrated using Seurat's own methods
# integration method, so we must remove it
print.status(paste("[INFO]: Reading in:", path.rds.1))
sobj.1 <- readRDS(path.rds.1)
sobj.1$region <- region.id.1 # extract from file name
DefaultAssay(sobj.1) <- "RNA"

sobj.1.integrated <- ("integrated" %in% names(sobj.1@assays))

if (sobj.1.integrated) {
    sobj.1[["integrated.old"]] <- sobj.1[["integrated"]]
    sobj.1[["integrated"]] <- NULL
}

print.status(paste("[INFO]: Reading in:", path.rds.2))
sobj.2 <- readRDS(path.rds.2)
sobj.2$region <- region.id.2 # extract from file name
DefaultAssay(sobj.2) <- "RNA"

sobj.2.integrated <- ("integrated" %in% names(sobj.2@assays))

if (sobj.2.integrated) {
    sobj.2[["integrated.old"]] <- sobj.2[["integrated"]]
    sobj.2[["integrated"]] <- NULL
}

# Visual inspection
p1 <- DimPlot(sobj.1, reduction = "umap", group.by="region", cols=c("#F8766D")) + 
        NoLegend() +
        ggtitle(region.id.1)
p2 <- DimPlot(sobj.2, reduction = "umap", group.by="region", cols=c("#00BFC4")) + 
        NoLegend() +
        ggtitle(region.id.2)
pg <- plot_grid(p1, p2)
p.list <- append(p.list, list(pg))


# Integrate ###################################################################
print.status("[INFO]: Integrating datasets.")
# setup before integration
# combine objects and preprocess
sobj.combined <- merge(x=sobj.1, y=sobj.2, add.cell.ids = c(region.id.1, region.id.2))

sobj.combined <- NormalizeData(sobj.combined) # normalizes the merged RNA assay
sobj.combined <- FindVariableFeatures(sobj.combined) # finds variable features in the merged RNA assay

sobj.combined <- RunFastMNN(object.list = Seurat::SplitObject(sobj.combined, split.by = 'orig.ident'))

# Run analysis workflow on integrated data
# N.B. We skip ScaleData and RunPCA, as this is not needed with RunFastMNN. That is because
# RunFastMNN produces another dim.reduction, "mnn", used for the workflow
sobj.combined <- RunUMAP(sobj.combined, reduction = "mnn", dims = 1:30)
sobj.combined <- FindNeighbors(sobj.combined, reduction = "mnn", dims = 1:30)
sobj.combined <- FindClusters(sobj.combined) # TODO: fancy resolution inference

# Visual inspection
p1 <- DimPlot(sobj.combined, pt.size=0.01, reduction = "umap", group.by="region") + ggtitle(paste(region.id.1, "and", region.id.2, "integrated - color by region"))
p2 <- DimPlot(sobj.combined, pt.size=0.01, reduction = "umap", label = TRUE) + ggtitle(paste(region.id.1, "and", region.id.2, "integrated - color by cluster"))
pg <- plot_grid(p1, p2)
p.list <- append(p.list, list(pg))

pg <- DimPlot(sobj.combined, reduction = "umap", split.by = "region", label = TRUE) + ggtitle(paste(region.id.1, "VS.", region.id.2))
p.list <- append(p.list, list(pg))


# MARKER DETECTION ############################################################
print.status("[INFO]: Running FindMarkers().")
# set identities and levels
Idents(sobj.combined) <- "seurat_clusters"
sobj.combined$celltype.region <- paste(Idents(sobj.combined), sobj.combined$region, sep='_')
sobj.combined$celltype <- Idents(sobj.combined)
Idents(sobj.combined) <- "celltype.region"

idents.all <- levels(sobj.combined)
idents.sobj.1 <- idents.all[grep(sobj.1@meta.data$region[1], idents.all)]
idents.sobj.2 <- idents.all[grep(sobj.2@meta.data$region[1], idents.all)]

# confirm levels
print(idents.all)
print(idents.sobj.1)
print(idents.sobj.2)

# region.1
sobj.1.markers <- FindMarkers(sobj.combined, ident.1=idents.sobj.1, ident.2=idents.sobj.2, only.pos=TRUE, verbose=FALSE)
sobj.1.markers$region <- region.id.1
sobj.1.markers$region_other <- region.id.2

# region.2
sobj.2.markers <- FindMarkers(sobj.combined, ident.1=idents.sobj.2, ident.2=idents.sobj.1, only.pos=TRUE, verbose=FALSE)
sobj.2.markers$region <- region.id.2
sobj.2.markers$region_other <- region.id.1

# Save results ################################################################
print.status("[INFO]: Saving output to disk.")

write.csv(sobj.1.markers, "d14D-nomux-markers_region.csv") # e.g. d14D-nomux-markers_regional.csv
write.csv(sobj.2.markers, "d14V-nomux-markers_region.csv")

saveRDS(sobj.combined, out.rds) # save combined data for plotting later

# save plots
path.plots.out <- paste0(dirname(out.rds), '/', integrated.id, "-plots.pdf")
raster_pdf(file=path.plots.out, 8, 7, res=200)
p.list
dev.off()

