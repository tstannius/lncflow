suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(stringr))



#' Filter out mitochondrial and ribosomal genes from seurat object
#'
#' First identify mitochondrial (^MT), ribosomal (^RPL and ^RPS) genes,
#' then remove the relevant rows from counts, data and meta features
#'
#' @param seurat.object The seurat object to apply filter to
#' @return seurat.object The filtered object
#'
#' @author Oliver K. Møller, Tobias O. Stannius
filter.mito.ribo <- function(seurat.object) {
  # extract mitochondrial and ribosomal genes
  mito.ribo.genes <- c(grep(pattern = "^MT-", x = rownames(x = seurat.object@assays$RNA@meta.features), value = T, ignore.case = T), 
                     grep(pattern = "^RPL", x = rownames(x = seurat.object@assays$RNA@meta.features), value = T, ignore.case = T),
                     grep(pattern = "^RPS", x = rownames(x = seurat.object@assays$RNA@meta.features), value = T, ignore.case = T))
  
  genes.to.use <- rownames(seurat.object@assays$RNA@meta.features)[!(rownames(seurat.object@assays$RNA@meta.features)
                                                                  %in% mito.ribo.genes)]

  # filter object
  seurat.object@assays$RNA@counts <- seurat.object@assays$RNA@counts[genes.to.use,]
  seurat.object@assays$RNA@data <- seurat.object@assays$RNA@data[genes.to.use,]
  seurat.object@assays$RNA@meta.features <- seurat.object@assays$RNA@meta.features[genes.to.use,]

  n.genes <- length(mito.ribo.genes)
  print(paste("filter.mito.ribo:", n.genes, "genes removed."))

  return(seurat.object)
}



#' Common Seurat workflow
#'
#' Perform scale, normalize, regress, PCA, UMAP, etc.
#'
#' @param seurat.object The seurat object to run workflow on
#' @return seurat.object Updated seurat object
#'
#' @author Tobias O. Stannius
common.workflow <- function(seurat.object) {
    stop("Not implemented yet")
    return(seurat.object)
}



#' Calculate silhouette index across all clusters in 1 resolution
#'
#' @param data : data.frame / vector containing cluster assignments
#'               for a given resolution
#' @param resolution : float, the resolution to compute Silhouette of
#'
#' @return mean silhouette width
#'
#' @author Martina Cardinali
silhouette_analysis <- function(data, resolution) {
    
    # calculate distances between cells -> cell.embeddings stores the coordinates for each cell in low-dimensional space
    distances <- dist(Embeddings(data, reduction = "pca"))
    # extract clusters, need to be integers
    clusters <- as.integer(as.character(resolution))
    # run silhouette analysis
    sil <- cluster::silhouette(clusters, dist = distances)
    # save average silhouette width
    avg_sil_width <- mean(sil[1:nrow(sil),3])
    return(avg_sil_width)
}



#' Calculate silhouette index for all resolutions in Seurat object
#'
#' @param sobj : Seurat Obect
#'
#' @return mean silhouette width
#'
#' @author Martina Cardinali, Tobias O. Stannius
silhouette_wrapper <- function(sobj) {
  cluster.assignments <- sobj@meta.data %>% select(matches("_snn_res.")) # extract resolutions
  resolutions <- str_remove(names(cluster.assignments), "(integrated|RNA|SCT)_snn_res.")

  # iterate over each resolution column containing cluster assignments and 
  # calculate the average silhouette width.
  # apply(X, MARGIN, FUN)
  # X is data, MARGIN=2 is to perform on columns, FUN is function
  all_avg_sil <- data.frame(avg_sil_width = apply(cluster.assignments, 2, function(x)
                      silhouette_analysis(sobj, x)), 
                      res = resolutions)
  return(all_avg_sil)
}



#' Plot silhouette index results
#'
#' @param data : data.frame
#' @param title : string
#'
#' @return mean silhouette width
#'
#' @author Martina Cardinali
plot.silhouette <- function(data, title){     # day_and_axis variable is used just for the plot title 
  ggplot(data=data, aes(x=res, y=avg_sil_width, group=1)) +
  geom_line(color="blue")+
  geom_point() +
  ggtitle(paste(title, "Silhouette")) +
  xlab("Resolution") +
  ylab("Average silhouette width") +
  theme_bw()
}



#' Query Seurat object for list of genes
#'
#' @param seurat.object : Seurat Object
#' @param gene.list     : list
#' @param threshold     : numeric
#'
#' @return list of expressed lncRNAs found in the Seurat Object
#'
#' @author Oliver K. Møller, Tobias O. Stannius
genes.in.sobj <- function(seurat.object, gene.list, threshold=100.0) {
  tmp.assay <- GetAssay(seurat.object)
  tmp.assay.sum <- Matrix::rowSums(tmp.assay) # sum over rows, i.e. genes
  tmp.assay.nonZero <- names(tmp.assay.sum[tmp.assay.sum > threshold])
  
  res <- intersect(gene.list, tmp.assay.nonZero)
  return(res)
}



#' Prefix msg with date-time stamp and print
#'
#' @param msg : str
#'
#' @return NULL
#'
#' @author Tobias O. Stannius
print.status <- function(msg) {
  print(paste(format(Sys.time(), "%y-%m-%d %H:%M:%S"), msg))
}
