# About

Pipeline for processing raw 10X scRNA-seq data to obtain `.loom` and `.rds` files, as well as running the default Seurat analysis pipeline.

### Pipeline

#### Generate count matrices
There are two tools available for generating count matrices: `Cellranger` by 10XGenomics and `Velocyto` by the Linnarsson Lab. The latter includes both intronic and exonic reads when generating count matrices, which is needed for counting lncRNAs.

##### Generate `rds`using Cellranger.
TODO

##### Generate `.loom` using Velocyto. Use `run_velocyto.sh` script.
Just provide experiment name - the script will look for the experiment name in `/nfsdata/data/data-runs/170907-kirkeby-mistr/`

Example
```
bash run_velocyto.sh d5d-5000_cells
```

##### Generate `.rds` from `.loom`. Use `loom2rds.R` script.

Example
```
./loom2rds.R /scratch/tstannius/velocyto/loom-files/d5d-5000_cells.loom
```

This will save a Seurat object in `out/d5d-5000_cells_raw`.


##### Generate HTO count matrice using CITE-seq
This is needed for the demultiplexing step.

##### Demultiplex and analyze data using Seurat. See template in `seurat_pipeline`.


##### Downstream analysis
E.g. find markers, dimensionality reduction or integration of datasets.

# Datasets

##### Day 14 dorsal

This was the first experiment to be sequenced and it should be treated differently for various reasons.

1. The day 14 MISTR tissue was supposed to be split in 5 parts, A-E, and sequenced separately, since we were not aware of HTO's. Thus it is not necessary to run CITE-seq count and do demultiplexing. (QUA: What about other sorts of QC?)
2. The day 14 MISTR tissue subregions were accidentally mixed. Originally, there should have been 5 runs (A-E), but A and B were pooled with C. Thus we have c, d and e.

The suggested strategy for the day 14 datasets is instead to integrate them using the Seurat approach (QUA: Use RNA or integrated assay?). 

