
# dev notes
# Find a generalized way to generate the time plot so that it's not hard coded
# 
# Things to try:
# Use ggrastr and see if it's faster: https://github.com/VPetukhov/ggrastr
# Plot from H5Seurat files and see if it's faster. Perhaps it's possible to make this with multiprocessing?

library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)
library(rasterpdf)
library(cowplot)



lnc_rna_list <- scan("/nfsdata/projects/tstannius/data/resources/hg38/ensembl/extracted/genes_lncRNA.txt", character(), quote="")

dir.data <- "/nfsdata/projects/tstannius/lncflow/out/hg38"
files <- list.files(dir.data, pattern = "(*seurat_preproc.rds)|(*markers.csv)", recursive = TRUE, full.names = TRUE)
files <- str_sort(files, numeric=TRUE)

filter.regions <- c("es", "d", "v")

experiments <- list()
markers.list <- list()


# iterate over files and store in nested list structure
c <- 1
for (f in files) {
  # extract experiment conditions from filename, e.g. day 14 ventral
  n <- tools::file_path_sans_ext(basename(f))
  run.id <- str_extract(n, pattern="[0-9]+(es|[dD]|[vV])(-|_)") # use - to ensure we don't get e.g. d14Dc
  run.id <- str_extract(run.id, pattern="[0-9]+(es|[dD]|[vV])")
  timepoint <- str_extract(run.id, pattern="[0-9]+")
  region <- tolower(str_extract(run.id, pattern="(es|[dD]|[vV])"))
  
  if (is.na((run.id))) {
    next
  }
  
  # determine if f is a rds or csv file
  key.file <- ""
  if (grepl("seurat_preproc.rds", f, fixed = TRUE)) {
    key.file <- "sobj"
    value <- readRDS(f)
    DefaultAssay(value) <- "RNA"
    if (run.id == "d14D") {
      DefaultAssay(value) <- "integrated"
    }
  } else if (grepl("markers.csv", f, fixed = TRUE)) {
    key.file <- "markers"
    markers <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
    markers$run_id <- run.id
    markers.list[[c]] <- markers
    value <- f
  } else {
    print("Unknown file type.")
    stop()
  }

  # if timepoint not already in experiment list, add it
  if (!timepoint %in% names(experiments)) {
    experiments[[timepoint]] <- list()
  }

  # if region not already in experiment list, add it
  if (!region %in% names(experiments[[timepoint]])) {
    experiments[[timepoint]][[region]] <- list()
  }

  # finally save file path
  experiments[[timepoint]][[region]][[key.file]] <- value

  c <- c + 1
}

# time.list <- str_sort(names(experiments), numeric = TRUE)

markers.combined <- dplyr::bind_rows(markers.list)
markers.combined <- markers.combined[row.names(markers.combined) %in% lnc_rna_list, ] # subset to lncRNA



feature_element <- function(run.id, sobj, feature) {
  # if expression is 0 in all cells, Seurat will throw a warning and color all cells purple
  # to prevent that, we force it to color all gray, if a warning is thrown
  p <- tryCatch({
    p <- FeaturePlot(sobj, feature, pt.size=0.01)
  }, warning = function(w) {
    p <- FeaturePlot(sobj, feature, cols = c("gray","gray"), pt.size=0.01)
  })
  
  p <- p + 
    NoLegend() + 
    ggtitle(run.id) + 
    NoAxes() # + small_font_theme
  return(p)
}

# feature_element("14v", experiments$`14`$v$sobj, "AC003975.1") # all should be gray
# feature_element("14v", experiments$`14`$v$sobj, "MALAT1") # all should be blue

TimePlotFeatures <- function(exp.list, feature.run.id, feature) {
  time.list <- str_sort(names(exp.list), numeric = TRUE)
  # time.list <- c("0", "1", "2", "9", "14", "35") # TODO: remove this part
  n.timepoints <- length(time.list)

  
  placeholder <- textGrob('')
  # handle special case 0, if present
  p0 <- NULL
  if ("0" %in% time.list) {
    time.list <- time.list[2:length(time.list)]
    sobj <- exp.list[["0"]][["es"]]$sobj
    run.id <- "0es"
    p0 <- feature_element(run.id, sobj, feature)
    p0 <- arrangeGrob(placeholder, p0, placeholder, ncol = 1, heights=c(.25,.5,.25))
  }
  
  ps <- lapply(time.list, FUN = function(t) {
    # dorsal
    sobj <- exp.list[[t]][["d"]]$sobj
    run.id <- paste0(t, "d")
    p.d <- feature_element(run.id, sobj, feature)
    
    # ventral
    sobj <- exp.list[[t]][["v"]]$sobj
    run.id <- paste0(t, "v")
    p.v <- feature_element(run.id, sobj, feature)
    
    return(arrangeGrob(p.d, p.v, ncol=1))
  })
  
  mask <- markers.combined$run_id==feature.run.id & markers.combined$gene==feature
  marker.meta <- markers.combined[mask,]
  # title <- paste(feature, "- marker in run", marker.meta$run_id, "cluster:", marker.meta$cluster)
  title <- feature
  ps <- c(list(p0), ps)
  ps <- arrangeGrob(grobs=ps, ncol=n.timepoints, top=title)
  
  return(ps)
}

# grid.arrange() actually shows the plot
# grid.arrange(TimePlotFeatures(experiments, "35D", "CDH6"))



TimePlotViolins <- function(exp.list, feature) {
  time.list <- str_sort(names(exp.list), numeric = TRUE)
  
  df <- data.frame()
  
  for (t in time.list) {
    regions <- names(experiments[[t]])

    for (r in regions) {
      condition <- paste0(t, r) # e.g. 1d

      tryCatch({
        df_query <- Seurat::FetchData(exp.list[[t]][[r]]$sobj, feature)
        df_query <- cbind(day = condition, df_query)
        colnames(df_query) <- c("day", "gene") # ggplot does not play well with odd chars in the name, so this is safer
  
        df <- rbind(df, df_query)
      },
      error = function(e){
        print(paste0("Error in dataset day: ", d, " ", e))
      })
    }
    
  }
  
  # make day an ordered factor
  df$day <- factor(df$day, levels=str_sort(unique(df$day), numeric = TRUE))
  p <- ggplot(df, aes(x=day, y=gene, fill=factor(day))) +
    geom_violin(trim=TRUE, scale="width") + # scale can be set to "area", "count" or "width"
    # geom_point(position="jitter") +
    labs(x="Condition", y="Expression") +
    theme_classic() +
    theme(legend.position="none") +
    ggtitle(feature)

  return(p)
}

# TimePlotViolins(experiments, "CDH6") # or LINC01405 or MALAT1



TimePlotDetails <- function(exp.list, run.id, feature) {
  timepoint <- str_extract(run.id, pattern="[0-9]+")
  region <- tolower(str_extract(run.id, pattern="(es)|[dD]|[vV]"))

  p1 <- FeaturePlot(exp.list[[timepoint]][[region]]$sobj, c(feature))

  p2 <- DimPlot(exp.list[[timepoint]][[region]]$sobj, label = TRUE) + ggtitle('')
  
  p.feature <- arrangeGrob(p1, p2, ncol=2)
}

# grid.arrange(TimePlotDetails(experiments, "35D", "CDH6"))


today <- format(Sys.time(), "%y%m%d")
path.plots <- paste0("./timeline_plots-lncRNA-", today, ".pdf")
n.genes <- dim(markers.combined)[1]
# raster_pdf(file=path.plots, width=10, height=15, res=200)
# for (i in 1:n.genes) {
#   g <- markers.combined$gene[i]
#   rid <- markers.combined$run_id[i]
#   print(paste(i, '/', n.genes, '-', rid, g))
# 
#   p.time.features <- TimePlotFeatures(experiments, rid, g)
#   p.time.violins <- TimePlotViolins(experiments, g)
#   p.feature <- TimePlotDetails(experiments, rid, g)
#   grid.arrange(p.time.features, p.time.violins, p.feature, nrow=3)
# }
# dev.off()

g <- "SIX3OS"
rid <- "999"
fp1pdf <- paste0("~/plots/timeline-lncRNA-feat-", g, ".pdf")
fp1png <- paste0("~/plots/timeline-lncRNA-feat-", g, ".png")
p.time.features <- TimePlotFeatures(experiments, rid, g)
p <- grid.arrange(p.time.features)
pdf(file=fp1pdf, width=10, height=5)
grid.arrange(p.time.features)
dev.off()
ggsave2(fp1png, p, width=10, height=5)

fp1pdf <- paste0("~/plots/timeline-lncRNA-vln-", g, ".pdf")
fp1png <- paste0("~/plots/timeline-lncRNA-vln-", g, ".png")
p.time.violins <- TimePlotViolins(experiments, g)
p <- grid.arrange(p.time.violins)
pdf(file=fp1, width=10, height=4)
grid.arrange(p.time.violins)
dev.off()
ggsave2(fp1png, p, width=10, height=4)

