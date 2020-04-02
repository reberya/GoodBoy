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
#                           - STAR
#                           - subread
#
# RUN: snakemake -s code/myworkflow.smk --cores 1 --latency-wait 30
################################################################################

from os.path import join
import os
import glob


##############################################################
# User Modifiables
##############################################################

# working directory
DATA_DIR = "/home/rebernrj/testenv/"

##############################################################
# Globals
##############################################################

# samples/read names
READS = ["1", "2"]
SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob(DATA_DIR + 'data/FQ/*.R1.fastq.gz')]

# directories
GENOME_DIR = DATA_DIR + "genome/GRCm38_vM21/"
FQC_DIR = DATA_DIR + "output/fastqc/"
STAR_DIR = DATA_DIR + "output/star/"
FEATURECOUNTS_DIR = DATA_DIR + "output/featureCounts/"


##############################################################
# List of end point files needed
##############################################################

FQC = expand(DATA_DIR + "output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = [DATA_DIR + "genome/GRCm38_vM21/genomeParameters.txt"]
ALIGNED = expand(DATA_DIR + "output/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
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
        files = expand(DATA_DIR + "data/FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS)
    output: FQC
    params: time = "01:00:00", mem = "1000m",
        fqc = FQC_DIR
    threads: 2
    shell: """ \
    mkdir -p {params.fqc}; \
    fastqc -o {params.fqc} {input.files} \
    """


# STAR genome index
rule starGenomeIndex:
    input:
        gd = GENOME_DIR
    output:
        index = GENOME_DIR + "/genomeParameters.txt"
    params: time="10:00:00", mem ="8000m", readLength="51",
        genome = GENOME_DIR
    threads: 6
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.gd} \
        --genomeFastaFiles {params.genome}GRCm38.primary_assembly.genome.fa \
        --sjdbGTFfile {params.genome}gencode.vM21.annotation.gtf \
        --sjdbOverhang {params.readLength}; \
        mv Log.out log/hpc/
        """


#STAR aligner
rule star:
    input:
        files = expand(DATA_DIR + "data/FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS),
        index = GENOME_DIR + "genomeParameters.txt"
    output: ALIGNED
    params: time="10:00:00", mem ="6000m",
        ext = expand(DATA_DIR + "output/star/{sample}_", sample=SAMPLES),
        gd = GENOME_DIR,
        star = STAR_DIR
    threads: 6
    shell:
        """
        mkdir -p {params.star}; \
        STAR \
        --genomeDir {params.gd} \
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
        bam = ALIGNED,
        gtf = GENOME_DIR + "gencode.vM21.annotation.gtf"
    output: FC
    params: time="10:00:00", mem="6000m",
        fc = FEATURECOUNTS_DIR,
        star = STAR_DIR
    threads: 2
    shell:
        """
        featureCounts -p -t exon -g gene_id \
        -a {input.gtf} \
        -o {params.fc}/counts.txt \
        {params.star}*Aligned.sortedByCoord.out.bam
        """
