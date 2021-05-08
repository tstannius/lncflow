#!/bin/bash

# usage: bash make_intronic-exonic_dataset.sh

# n.b. this script produces a loom file. To use this file with R and Seurat, use the loom_to_rds.R script

# todo turn this into a snakemake pipeline that supports creating the outputs that dont exist from the input

set -euo pipefail

# declare -a EXPERIMENT_NAMES=() # this must be set, e.g. to "d5d" # declare array not necessary for this simplified version
# furthermore, this would run everything sequentially, which is not ideal
EXPERIMENT_NAME=$1 # e.g. "d5d-5000_cells"
echo "Running: ${EXPERIMENT_NAME}"

#module use  /usr/local/mtools/modules/
#module load anaconda/4.8.2
source activate velocyto

# constants
DATARUNS_DIR=/nfsdata/data/data-runs
REFDATA_DIR=/scratch/tstannius/RefData # should be in a common location
PROJECT_ID='170907-kirkeby-mistr'
ANNOTATION_FILE='Homo_sapiens.GRCh37.82.gtf'
REPEAT_MASKER='hg19_rmsk.gtf'
WORKING_DIR=/scratch/${USER}/velocyto

# make workdir
mkdir -p ${WORKING_DIR}

# copy 10X directory to workdir on scratch in order to have writing permission
cp -r ${DATARUNS_DIR}/${PROJECT_ID}/${EXPERIMENT_NAME} ${WORKING_DIR}/
chmod g+w ${WORKING_DIR}/${EXPERIMENT_NAME}

# run velocyto
# Q: how does this become intronic + exonic?
# A: velocyto is designed to also recognize and pick up unspliced mRNA
# Usage: velocyto run10x [OPTIONS] SAMPLEFOLDER GTFFILE
velocyto run10x \
    --mask ${REFDATA_DIR}/${REPEAT_MASKER} \
    --samtools-threads 40 \
    ${WORKING_DIR}/${EXPERIMENT_NAME} \
    ${REFDATA_DIR}/${ANNOTATION_FILE} \
    -v \

# run velocyto with advanced settings
# n.b. the only difference here seems to be the use of the --bcfile
# what does the barodes file do? 
# velocyto run \
#     --bcfile ${WORKING_DIR}/${EXPERIMENT_NAME}/outs/filtered_gene_bc_matrices/hg19/barcodes.tsv \
#     --outputfolder ${WORKING_DIR}/${EXPERIMENT_NAME}/velocyto \
#     --sampleid ${EXPERIMENT_NAME} \
#     --mask ${REFDATA_DIR}/${REPEAT_MASKER} \
#     ${WORKING_DIR}/${EXPERIMENT_NAME}/outs/possorted_genome_bam.bam \
#     ${REFDATA_DIR}/${ANNOTATION_FILE}

# rescue .loom file
mkdir -p ${WORKING_DIR}/loom-files/
mv ${WORKING_DIR}/${EXPERIMENT_NAME}/velocyto/${EXPERIMENT_NAME}.loom ${WORKING_DIR}/loom-files/

# record session info
pip freeze > ${WORKING_DIR}/loom-files/${EXPERIMENT_NAME}_velocyto_envinfo.txt

# Free some space
rm -r ${WORKING_DIR}/${EXPERIMENT_NAME}
