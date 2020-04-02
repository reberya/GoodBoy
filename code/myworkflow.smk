################################################################################
# TITLE:                    2020 RNAseq pipeline
# AUTHOR:                   Ryan Rebernick
# DATE LAST MODIFIED:       03/30/2020
#
# FUNCTION:                 FastQC
#                           STAR alignment
#
#
# USES:
#                           - python3.7-anaconda/2019.07
#                           - smk_env1
#                           - fastQC
#
# RUN: snakemake -s code/myworkflow.smk --cores 1 --latency-wait 30
################################################################################

from os.path import join
import os
import glob
import numpy as np


##############################################################
# Globals
##############################################################

# samples/read names
READS = ["1", "2"]
SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob('/home/rebernrj/testenv/data/FQ/*.R1.fastq.gz')]
GENOMEDIR = "/home/rebernrj/testenv/genome/GRCm38_vM21/"

##############################################################
# List of directories needed and end point files for analysis
##############################################################

FQC = expand("/home/rebernrj/testenv/output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = ["/home/rebernrj/testenv/genome/GRCm38_vM21/genomeParameters.txt"]
ALIGNED = expand("/home/rebernrj/testenv/output/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)


##############################################################
# Rules
##############################################################


# end files required
rule all:
    input: FQC + GENOME + ALIGNED
    params: time="10:00:00", mem="50m"


# FastQC
rule fastqc:
    input:
        files = expand("/home/rebernrj/testenv/data/FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS)
    output:
        reads = expand("/home/rebernrj/testenv/output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
    params: time="01:00:00", mem="1000m"
    threads: 2
    shell: """ \
    mkdir -p /home/rebernrj/testenv/output/fastqc; \
    fastqc -o /home/rebernrj/testenv/output/fastqc/ {input.files} \
    """

# STAR genome index
rule starGenomeIndex:
    input:
        gd = GENOMEDIR
    output:
        index = "/home/rebernrj/testenv/genome/GRCm38_vM21/genomeParameters.txt"
    params: time="10:00:00", mem="8000m", readLength="51"
    threads: 6
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.gd} \
        --genomeFastaFiles /home/rebernrj/testenv/genome/GRCm38_vM21/GRCm38.primary_assembly.genome.fa \
        --sjdbGTFfile /home/rebernrj/testenv/genome/GRCm38_vM21/gencode.vM21.annotation.gtf \
        --sjdbOverhang {params.readLength}; \
        mv Log.out log/hpc/
        """

#STAR aligner
rule star:
    input:
        files = expand("/home/rebernrj/testenv/data/FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS),
        gd = GENOMEDIR,
        index = "/home/rebernrj/testenv/genome/GRCm38_vM21/genomeParameters.txt"
    output:
        bam = expand("/home/rebernrj/testenv/output/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    params: time="10:00:00", mem="6000m", ext = expand("/home/rebernrj/testenv/output/star/{sample}/", sample=SAMPLES)
    threads: 6
    shell:
        """
        mkdir -p /home/rebernrj/testenv/output/star; \
        STAR \
        --genomeDir {input.gd} \
        --runThreadN {threads} \
        --readFilesIn {input.files} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.ext} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard
        """
