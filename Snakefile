#!/usr/local/mtools/tools/anaconda/4.8.2/bin/snakemake -s

import os
import builtins
import re
from pathlib import Path

# vars  #######################################################################
# TODO: migrate to config file
env_velocyto = "/nfsdata/projects/tstannius/data/resources/envs/velocyto/0.17"
env_r = "/nfsdata/projects/tstannius/data/resources/envs/R/4.0.3"
hgref_version = "hg38"
path_annotation = f"/nfsdata/projects/tstannius/data/resources/hg38/annotation.gtf"
path_rmsk = f"/nfsdata/projects/tstannius/data/resources/hg38/rmsk.gtf"
path_cellranger_out = "/nfsdata/data/data-runs/mistr_cellranger5.0/170907-kirkeby-mistr"
path_citeseq_out = "/nfsdata/projects/tstannius/data/citeseq"
path_workdir = "/nfsdata/projects/tstannius/lncflow"
path_scratch = "/scratch/tstannius/lncflow"


# set up env ##################################################################
workdir: path_workdir

# snakemake normally makes dirs, but because we want to use symlinks, they must
# be set up first. We need the symlink setup so that we can keep loom files 
# on scratch, but pretend they are local
os.makedirs(f"{path_workdir}/out/{hgref_version}", exist_ok=True)
os.makedirs(f"{path_scratch}/out/{hgref_version}/loom", exist_ok=True)
os.makedirs(f"{path_scratch}/tmp", exist_ok=True)
s1 = Path(f"{path_scratch}/out/{hgref_version}/loom")
# s2 = Path(f"{path_scratch}/tmp")
if not s1.exists():
    os.symlink(src=s1, dst=f"{path_workdir}/out/{hgref_version}/loom")
# if not s2.exists():
#     os.symlink(src=s2, dst=f"{path_workdir}/tmp")

# define targets
# runs - TODO: migrate to sample sheet
runs = [
        "d0es",
        "d1d",
        "d1v",
        "d2D",
        "d2V",
        "d5D",
        "d5V",
        # "d9D", # bad run
        # "d9V", # bad run
        "d9D_2", # fix poor naming
        "d9V_2",
        "d14V", 
        "d14D",
        "d14Dc", 
        "d14Dd", 
        "d14De",
        "d35D",
        "d35V"
        ]

days = set([re.search("\d+", r).group(0) for r in runs if r != "d0es"])

runs_citeseq = [name for name in os.listdir(path_citeseq_out)]

runs2demux = ["demux" if r in runs_citeseq else "nomux" for r in runs]


workflow.first_rule = "all"

wildcard_constraints:
    run = r"|".join(set(runs)),
    mux = r"|".join(set(runs2demux)),
    day = r"|".join(set(days))


# input functions #############################################################
def get_loom(w, p):
    return f"{p}/tmp/{w}-5000_cells/velocyto/{w}-5000_cells.loom"

def get_cite_seq_dir(w, p):
    return f"{p}/{w}/umi_count"

def get_rds_region(w, region, hgref, runs, runs2demux):
    key = f"d{w}{region}"
    sample = None
    mux = None
    for r,m in zip(runs, runs2demux):
        if key.lower() in r.lower():
            val = r
            mux = m
            break
    return f"out/{hgref}/{val}/{val}-{mux}-seurat_preproc.rds"



# rules #######################################################################
rule copy_10x_samplefolder:
    input:
        f"{path_cellranger_out}/{{run}}-5000_cells"
    output:
        temp(directory(f"{path_scratch}/tmp/{{run}}-5000_cells")) # on scratch
    threads: 1
    shell:
        """
        cp -r {input} {output}
        """
# try sorting to fix issue with velocyto
# samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam
# alternatively do not use velocyto, as it is not throughly investigated if there's any benefit
# over using cellranger

rule velocyto:
    input:
        samplefolder = rules.copy_10x_samplefolder.output # TODO does this work? f"tmp/{{run}}-5000_cells"
        # samplefolder = "/scratch/tstannius/bicropflow/tmp/TEST_bicrop_mRNA" # For gRNA detection pilot
    output:
        f"out/{hgref_version}/loom/{{run}}.loom" # on scratch
        # f"out/{hgref_version}/loom/d0-dv-pgRNA.loom" # on scratch
    params:
        loom = lambda wildcards: get_loom(wildcards, path_scratch),
        # loom = f"{path_scratch}/tmp/TEST_bicrop_mRNA/velocyto/TEST_bicrop_mRNA.loom", # For gRNA detection pilot
        conda_env = env_velocyto,
        gtffile = path_annotation,
        # gtffile = "/nfsdata/projects/tstannius/data/resources/hg38/cropseq/GRCh38-2020-A-spiked/genes/genes.gtf", # For gRNA detection pilot
        rmsk = path_rmsk,
        mem = 1500
    threads: 20
    log: f"logs/velocyto/{hgref_version}/{{run}}.log"
    # log: f"logs/velocyto/{hgref_version}/d0-dv-pgRNA.log" # For gRNA detection pilot
    shell:
        "("
        "module load anaconda;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "velocyto run10x"
        "    --mask {params.rmsk}"
        "    --samtools-threads {threads}"
        "    --samtools-memory {params.mem}"
        "    {input.samplefolder}"
        "    {params.gtffile};\n"

        # rescue loom
        "mv {params.loom} {output};\n"

        ") 2> {log}"



rule loom2rds:
    input:
        rules.velocyto.output
    output:
        temp(f"out/{hgref_version}/{{run}}/{{run}}-raw.rds") # on nfsdata
        # f"out/{hgref_version}/d0-dv-pgRNA/d0-dv-pgRNA-raw.rds" # for gRNA detection pilot
    log: f"logs/loom2rds/{hgref_version}/{{run}}.log"
    # log: f"logs/loom2rds/{hgref_version}/d0-dv-pgRNA-raw.log" # for gRNA detection pilot
    params:
        conda_env = env_r,
        path_script = "/nfsdata/projects/tstannius/lncflow/scripts/seurat/loom2rds_seurat4.R" # abs path necessary
    threads: 4
    shell:
        "("
        "module load anaconda hdf5;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "Rscript {params.path_script} -i {input} -o {output};\n"

        ") 2> {log}"



rule demux:
    input:
        rds = rules.loom2rds.output,
        cite_seq_dir = lambda wildcards: get_cite_seq_dir(wildcards, path_citeseq_out)
    output:
        f"out/{hgref_version}/{{run}}/{{run}}-demux.rds"
    log: f"logs/demux/{hgref_version}/{{run}}.log"
    params:
        conda_env = env_r,
        path_script = "/nfsdata/projects/tstannius/lncflow/scripts/seurat/demultiplex.R" # abs path necessary
    threads: 4
    shell:
        "("
        "module load anaconda hdf5;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "Rscript {params.path_script} -r {input.rds} -c {input.cite_seq_dir} -o {output};\n"

        ") 2> {log}"



rule nomux:
    """Dummy rule for samples without demultiplexing
    """
    input:
        rules.loom2rds.output
    output:
        f"out/{hgref_version}/{{run}}/{{run}}-nomux.rds"
    threads: 1
    shell:
        "cp {input} {output}"



rule seurat_preprocess:
    input:
        f"out/{hgref_version}/{{run}}/{{run}}-{{mux}}.rds"
    output:
        f"out/{hgref_version}/{{run}}/{{run}}-{{mux}}-seurat_preproc.rds"
    log: f"logs/seurat_preprocess/{hgref_version}/{{run}}-{{mux}}.log"
    params:
        conda_env = env_r,
        path_script = "/nfsdata/projects/tstannius/lncflow/scripts/seurat/preprocess.R", # abs path necessary
        path_dirout = f"out/{hgref_version}/{{run}}"
    threads: 4
    shell:
        "("
        "module load anaconda;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "Rscript {params.path_script} -i {input} -o {params.path_dirout} -l;\n"

        ") 2> {log}"



rule seurat_find_markers:
    input:
        rules.seurat_preprocess.output
    output:
        f"out/{hgref_version}/{{run}}/{{run}}-{{mux}}-markers.csv"
    log: f"logs/seurat_find_markers/{hgref_version}/{{run}}-{{mux}}.log"
    params:
        conda_env = env_r,
        path_script = "/nfsdata/projects/tstannius/lncflow/scripts/seurat/find_markers.R", # abs path necessary
    threads: 4
    shell:
        "("
        "module load anaconda;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "Rscript {params.path_script} -i {input} -o {output};\n"

        ") 2> {log}"



rule seurat_find_markers_regional:
    input:
        rds1 = lambda wildcards: get_rds_region(wildcards, 'D', hgref_version, runs, runs2demux), # f"out/{hgref_version}/d{{day}}D/d{{day}}-seurat_preproc.rds",
        rds2 = lambda wildcards: get_rds_region(wildcards, 'V', hgref_version, runs, runs2demux) # f"out/{hgref_version}/d{{day}}V/d{{day}}-seurat_preproc.rds"
    output:
        csv1 = f"out/{hgref_version}/d{{day}}/d{{day}}D-markers_regional.csv",
        csv2 = f"out/{hgref_version}/d{{day}}/d{{day}}V-markers_regional.csv",
        rds = f"out/{hgref_version}/d{{day}}/d{{day}}-integrated.rds"
    log: f"logs/seurat_find_markers_regional/{hgref_version}/d{{day}}-find_markers_regional.log"
    params:
        conda_env = env_r,
        path_script = "/nfsdata/projects/tstannius/lncflow/scripts/seurat/find_markers_regional.R", # abs path necessary
    threads: 4
    shell:
        "("
        "module load anaconda;\n"
        "set +u; source activate {params.conda_env}; set -u;\n"

        "Rscript {params.path_script}"
        "   --input_rds_1 {input.rds1}"
        "   --input_rds_2 {input.rds2}"
        "   --output_csv_1 {output.csv1}"
        "   --output_csv_2 {output.csv2}"
        "   --output_rds_integrated {output.rds};\n"

        ") 2> {log}"


rule all:
    input:
        expand("out/{hgrv}/{sample}/{sample}-{mux}-markers.csv", zip, hgrv=[hgref_version for i in range(len(runs))], sample=runs, mux=runs2demux),
        # expand("out/{hgrv}/d{day}/d{day}{region}-markers_regional.csv", hgrv=[hgref_version for i in range(len(days))], day=days, region=["D", "V"]),
        # expand("out/{hgrv}/d{day}/d{day}-integrated.rds", hgrv=[hgref_version for i in range(len(days))], day=days),
        


