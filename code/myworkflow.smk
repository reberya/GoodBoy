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
GENOME_DIR = "/home/rebernrj/testenv/genome/GRCm38_vM21/"
STAR_DIR = "/home/rebernrj/testenv/output/star/"
FEATURECOUNTS_DIR = "/home/rebernrj/testenv/output/FeatureCounts/"

##############################################################
# List of directories needed and end point files for analysis
##############################################################

FQC = expand("/home/rebernrj/testenv/output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = ["/home/rebernrj/testenv/genome/GRCm38_vM21/genomeParameters.txt"]
ALIGNED = expand("/home/rebernrj/testenv/output/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
FC = [FEATURECOUNTS_DIR + "counts.txt"]

##############################################################
# Rules
##############################################################


# end files required
rule all:
    input: FQC + GENOME + ALIGNED + FC
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
        gd = GENOME_DIR
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
        gd = GENOME_DIR,
        index = "/home/rebernrj/testenv/genome/GRCm38_vM21/genomeParameters.txt"
    output:
        bam = expand("/home/rebernrj/testenv/output/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        dir = STAR_DIR
    params: time="10:00:00", mem="6000m", ext = expand("/home/rebernrj/testenv/output/star/{sample}_", sample=SAMPLES)
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

# Feature counts
rule featureCounts:
    input:
        dir = STAR_DIR,
        gtf = GENOME_DIR + "*.gtf"
    output: FC
    params: time="10:00:00", mem="6000m", dir = FEATURECOUNTS_DIR
    threads:
    shell:
    """
    mkdir -p /home/rebernrj/testenv/output/featureCounts; \
    featureCounts -p -t exon -g gene_id \
    -a {input.gtf} \
    -o {output.dir}/counts.txt \
    {input.dir}*Aligned.sortedByCoord.out.bam
    """
