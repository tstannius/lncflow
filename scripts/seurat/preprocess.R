#!/usr/bin/env Rscript

# usage
#   this script runs a default Seurat workflow
#   incl. preprocessing, QC, reports, UMAP, etc.

# based on the following vignettes and guides:
#   guided clustering:      https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#   cell cycle regression:  https://satijalab.org/seurat/articles/cell_cycle_vignette.html
#   SCTransform:            https://satijalab.org/seurat/articles/sctransform_vignette.html
#   qc tips:                https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html

# notes
# - none atm

# PARSE ARGUMENTS #############################################################
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))

option.list <- list(
    make_option(c("-i", "--input_rds"), type="character",
                  help="Input Seurat .rds file."),
    make_option(c("-o", "--outdir"), type="character",
                  help="Output directory."),
    make_option(c("-l", "--log"), action="store_true", default=FALSE,
                  help="Enable logging.")
)

args <- parse_args(OptionParser(option_list=option.list))

if (is.null(args$input_rds)) {
  stop("[ERROR]: No input .rds file provided. See `--help`.")
}
if (is.null(args$outdir)) {
  stop("[ERROR]: No output directory path provided. See `--help`.")
}

path.rds.in <- args$input
path.dir.out <- args$outdir
if (!str_ends(path.dir.out, '/')) {
    path.dir.out <- paste0(path.dir.out, '/')
}

project.id <- tools::file_path_sans_ext(basename(path.rds.in))

# list for storing qc plots
p.qc.list <- list()


# RUN WORKFLOW ################################################################
if (args$log) {
  # start log, if path provided
  path.log <- paste0(path.dir.out, project.id, "-preprocess.log")
  sink(path.log)
}

source("/nfsdata/projects/tstannius/lncflow/scripts/seurat/helpers.R")

print.status("[INFO]: Preprocessing initialized.")
print.status(paste("[INFO]: Input:", path.rds.in))
print.status("[INFO]: Loading dependencies.")
# if optparse completes, continue with rest of script
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(ggplot2))   # for qc report plots
suppressPackageStartupMessages(library(rasterpdf)) # for qc report plots
set.seed(42)

print.status("[INFO]: Loading Seurat rds.")
sobj <- readRDS(path.rds.in)

# set identity and assay
print.status("[INFO]: Setting Idents and DefaultAssay.")
sobj.integrated <- ("integrated" %in% names(sobj@assays))
if (!sobj.integrated) {
  Idents(sobj) <- "orig.ident"
  DefaultAssay(sobj) <- "RNA"
} else {
  print.status("[INFO]: Integrated dataset. Assuming dataset preprocssed and skipping QC.")
  DefaultAssay(sobj) <- "integrated"
}

print(sobj) # print info

# to keep track of how many removed in QC
n.cells.raw <- dim(sobj@meta.data)[1]

# QC and cell subsetting ######################################################
# TODO: Handle this elegantly, instead of one big if expression
if (!sobj.integrated) {

  print.status("[INFO]: Running QC and filtering.")
  # Common QC metrics
  # - The number of unique genes detected in each cell.
  #   - Low-quality cells or empty droplets will often have very few genes
  #   - Cell doublets or multiplets may exhibit an aberrantly high gene count
  # - Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
  # - The percentage of reads that map to the mitochondrial genome
  #   - Low-quality / dying cells often exhibit extensive mitochondrial contamination
  #   - We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
  #   - We use the set of all genes starting with MT- as a set of mitochondrial genes



  # Determine threshold for Mitochondrial gene content cutoff ###################
  print.status("[INFO]: Calculating mitochondrial gene content.")
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

  p.qc.vln.all <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  # determine best threshold for percent mt cutoff and make plot
  threshold.mt.q <- quantile(sobj$percent.mt, 0.90)
  threshold.mt.max <- 10
  threshold.mt.min <- 5
  threshold.mt.selected <- min(max(threshold.mt.min,threshold.mt.q),threshold.mt.max)
  print.status(paste("[INFO]: Selected threshold for percent.mt: ", threshold.mt.selected))

  # Visualize QC metrics as a violin plot
  p.qc.vln.percent_mt <- VlnPlot(sobj, features = c("percent.mt")) +
          geom_hline(yintercept=threshold.mt.max, color="red") +
          geom_hline(yintercept=threshold.mt.min, color="blue") +
          geom_hline(yintercept=threshold.mt.q, color="green")



  # Determine thresholds for low and too high feature counts ####################
  print.status("[INFO]: Calculating thresholds for nFeature_RNA.")
  threshold.feat.up <- quantile(sobj$nFeature_RNA, 0.99) # TODO: Increase to 0.95
  threshold.feat.lo <- quantile(sobj$nFeature_RNA, 0.01) # TODO: Increase to 0.05 or find a better way to handle d5V
  print.status(paste("[INFO]: 0.99 quantile threshold for nFeature_RNA:", threshold.feat.up))
  print.status(paste("[INFO]: 0.01 quantile threshold for nFeature_RNA:", threshold.feat.lo))

  # Visualize QC metrics as a violin plot
  p.qc.vln.nfeature_rna <- VlnPlot(sobj, features = c("nFeature_RNA")) +
      geom_hline(yintercept=c(threshold.feat.up), color="red") +
      geom_hline(yintercept=c(threshold.feat.lo), color="blue")


  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  p.qc.feat_scatter.rna.vs.mt <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p.qc.feat_scatter.count.vs.feat <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  

  # Subset based on thresholds ##################################################
  print.status("[INFO]: Subsetting dataset based on percent.mt and nFeature_RNA thresholds.")
  # Remove cells 
  # - with low or high feature content. 
  # Low indicates empty droplet or dying cell. 
  # High indicates dublet. See: https://www.biostars.org/p/407036/
  sobj <- subset(sobj, subset = 
                nFeature_RNA > threshold.feat.lo & 
                nFeature_RNA < threshold.feat.up & 
                percent.mt < threshold.mt.selected)

  print(sobj)
  n.cells.proc <- dim(sobj@meta.data)[1] # keep track of how many removed
  n.cells.removed <- n.cells.raw - n.cells.proc
  print.status(paste("Seurat::subset:", n.cells.removed, "cells removed."))

  # Remove genes
  print.status("[INFO]: Subsetting dataset to remove mitochondrial and ribosomal genes.")
  # i.e. mitochondrial and ribosomal genes to to avoid influence on downstream analysis
  sobj <- filter.mito.ribo(sobj)

  print(sobj)


  # Make plots after subsetting
  p.qc.vln.all.post <- VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p.qc.feat_scatter.rna.vs.mt.post <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p.qc.feat_scatter.count.vs.feat.post <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


  # Cell-cycle scoring, Normalize, ScaleData, and FindVariableFeatures ##########
  print.status("[INFO]: Running CellCycleScoring.")
  # CellCycleScoring requires data to be normalized, see:
  # https://github.com/satijalab/seurat/issues/3692
  sobj <- NormalizeData(sobj)

  # perform cell cycle scoring
  cc.genes <- cc.genes.updated.2019 # use updated set
  sobj <- CellCycleScoring(sobj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)

  # run SCTransform
  # Replaces Normalize, Scale and FindVariableFeatures
  # n.b. not recommended to regress out mitochondrial and ribosomal genes
  # 
  # Running SCTransform may yield a warning message, but this can be safely
  # ignored, according to: https://github.com/ChristophH/sctransform/issues/25#issuecomment-494674451

  print.status("[INFO]: Running SCTransform with cell cycle regression.")
  sobj <- SCTransform(sobj, vars.to.regress = c('S.Score', 'G2M.Score'), verbose = FALSE)
  print(sobj) # n.b. new assay, SCT, added



  # save all qc plots in list
  p.qc.list <- append(p.qc.list, list(p.qc.vln.all, 
              p.qc.vln.percent_mt, 
              p.qc.vln.nfeature_rna,
              p.qc.feat_scatter.rna.vs.mt,
              p.qc.feat_scatter.count.vs.feat,
              p.qc.vln.all.post,
              p.qc.feat_scatter.rna.vs.mt.post,
              p.qc.feat_scatter.count.vs.feat.post))
}

# Dimensionality reduction ####################################################
print.status("[INFO]: Running dimensionality reduction.")
# from sctransform vignette
sobj <- RunPCA(sobj, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)



# pca qc plot - no matter, we will use first 30 PC dims
p.qc.elbowplot <- ElbowPlot(sobj)



# Clustering ##################################################################
print.status("[INFO]: Running clustering.")
sobj <- FindNeighbors(sobj, dims = 1:30, verbose = FALSE) # use many dims for these complex datasets


# Find clusters using different clustering resolutions, until upper limit for n clusters
# first drop any previously calculated neighbors or clustering, e.g. if sobj is integrated
sobj@meta.data <- sobj@meta.data %>% select(-matches("_snn_res.|seurat_clusters"))
target.n.clusters <- 24
res.not.converged <- TRUE
res <- 0.4

while (res.not.converged) {
    print.status(paste("[INFO]: Clustering with resolution: ", res))

    sobj <- FindClusters(sobj, resolution = res, verbose = FALSE) # Seurat devs recommend 0.4-1.2 for 3K cells
    
    if (nlevels(sobj@meta.data$seurat_clusters) >= target.n.clusters) {
        print.status(paste("[INFO]: Converged on", target.n.clusters, "clusters with resolution:", res))
        res.not.converged <- FALSE
    } else {
        res <- res + 0.1
    }
}

print.status("[INFO]: Calculating best clustering resolution based on average silhouette width.")
# Compute Average Silhouette Width for different resolutions to determine best resolution
resolution.silhouettes <- silhouette_wrapper(sobj)

# qc plot
p.qc.resolution.silhouette <- ggplot(data=resolution.silhouettes, aes(x=res, y=avg_sil_width, group=1)) +
  geom_line(color="blue")+
  geom_point() +
  ggtitle(paste("Silhouette results of ", project.id)) +
  xlab("Resolution") +
  ylab("Average silhouette width") +
  theme_bw()

# select resolution that maximises Average Silhouette Width
resolution.max.idx <- which.max(resolution.silhouettes$avg_sil_width) # find index of max average silhouette width
resolution.best <- as.numeric(resolution.silhouettes$res[resolution.max.idx]) # find corresponding resolution to max avg silhouette width
print.status(paste("[INFO]: Resolution based on highest average silhouette width:", resolution.best))

# drop calculated resolutions and cluster assignments
sobj@meta.data <- sobj@meta.data %>% select(-matches("_snn_res.|seurat_clusters"))

print.status("[INFO]: Clustering with optimal resolution.")
# find clusters for best resolution
sobj <- FindClusters(sobj, resolution = resolution.best, verbose = FALSE) # Seurat devs recommend 0.4-1.2 for 3K cells



# Visualization ###############################################################
# finally plot umap with cluster assignments
p.dimplot.umap <- DimPlot(sobj, label = TRUE)

p.qc.list <- append(p.qc.list, list(p.qc.elbowplot,
               p.qc.resolution.silhouette,
               p.dimplot.umap))

# Save results ################################################################
print.status("[INFO]: Saving output to disk.")
# save rds
path.rds.out <- paste0(path.dir.out, project.id, "-seurat_preproc.rds")
saveRDS(sobj, path.rds.out)



# save qc plots
path.plots_qc.out <- paste0(path.dir.out, project.id, "-qc_plots.pdf")
raster_pdf(file=path.plots_qc.out, 8, 7, res=200)
p.qc.list
dev.off()



# save marker plots
path.plots_markers.out <- paste0(path.dir.out, project.id, "-marker_plots.pdf")

# plot each marker as feature plot
Idents(sobj) <- "orig.ident"
DefaultAssay(sobj) <- "RNA"
raster_pdf(file=path.plots_markers.out, 8, 7, res=200)
markers.regional <- sort(read.csv("/nfsdata/projects/tstannius/lncflow/scripts/seurat/markers.csv", stringsAsFactors = FALSE)$gene)

# if a gene is not expressed, print a warning, but don't throw error
for (m in markers.regional) {
  tryCatch({
    FeaturePlot(sobj, features=c(m), combine=F, label=TRUE)
  }, error = function(e) {
    print.status(paste("[WARNING]:", conditionMessage(e)))
    })
}
dev.off() # stop writing to pdf


print.status("[INFO]: Preprocessing complete.")

if (args$log) {
  # close log, if path provided
  sink()
}
