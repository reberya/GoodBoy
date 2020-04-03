################################################################################
# TITLE:                    2020 RNAseq pipeline (Goodboy)
# AUTHOR:                   Ryan Rebernick
# DATE LAST MODIFIED:       04/03/2020
#
# FUNCTION:                 FastQC
#                           STAR alignment
#                           FeatureCounts (subread)
#                           MultiQC
#
# USES:
#                           - python3.7-anaconda/2019.07
#                           - config/python_env/smk_env1.yaml
#                           - fastQC
#                           - STAR
#                           - subread
#                           - multiqc
#                           - config/python_env/multiqc_env1.yaml
#
# RUN: snakemake -s code/myworkflow.smk --cores 1 --latency-wait 30
################################################################################

from os.path import join
import os
import glob


##############################################################
# User Modifiables
##############################################################

# working directories
DATA_DIR = "/home/rebernrj/testenv/"
GOODBOY_DIR = "/home/rebernrj/GoodBoy/"
GENOME_DIR = "/home/rebernrj/testenv/genome/GRCm38_vM24"


##############################################################
# Globals
##############################################################

# samples/read names
READS = ["1", "2"]
SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob(DATA_DIR + 'FQ/*.R1.fastq.gz')]

# FASTQ files
FQ = expand(DATA_DIR + "FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS)

# directories
FQC_DIR = DATA_DIR + "output/fastqc/"
STAR_DIR = DATA_DIR + "output/star/"
FEATURECOUNTS_DIR = DATA_DIR + "output/featureCounts/"
MULTIQC_DIR = DATA_DIR + "output/multiqc/"
NEW_LOG_DIR = DATA_DIR + "output/logs/"
CONFIG_DIR = GOODBOY_DIR + "config/python_config/"
LOG_DIR = GOODBOY_DIR + "log/hpc/"



##############################################################
# List of end point files needed
##############################################################

FQC = expand(DATA_DIR + "output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)
GENOME = [DATA_DIR + "genome/GRCm38_vM21/genomeParameters.txt"]
ALIGNED = expand(DATA_DIR + "output/star/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
FC = [FEATURECOUNTS_DIR + "counts.txt"]
MULTIQC = [MULTIQC_DIR + "multiqc_report.html"]


##############################################################
# Rules
##############################################################

# end files required
rule all:
    input: FQC + GENOME + ALIGNED + FC + MULTIQC
    params: time="10:00:00", mem="50m",
        logDir = LOG_DIR,
        newLogDir = NEW_LOG_DIR
    shell:
        """
        mkdir -p {params.newLogDir}
        mv {params.logDir}* {params.newLogDir}
        """


# FastQC
rule fastqc:
    input: FQ
    output: FQC
    params: time = "01:00:00", mem = "1000m",
        fqc = FQC_DIR
    threads: 2
    shell: """ \
    mkdir -p {params.fqc}; \
    fastqc -o {params.fqc} {input} \
    """


# STAR genome index
rule starGenomeIndex:
    input:
        gd = GENOME_DIR
    output:
        index = GENOME
    params: time="10:00:00", mem ="8000m", readLength="51",
        genome = GENOME_DIR
    threads: 6
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {input.gd} \
        --genomeFastaFiles {params.genome}*.fa \
        --sjdbGTFfile {params.genome}*.gtf \
        --sjdbOverhang {params.readLength}; \
        mv Log.out log/hpc/
        """


#STAR aligner
rule star:
    input:
        files = FQ,
        index = GENOME
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


# multiqc
rule multiqc:
    input: FQC + GENOME + ALIGNED + FC
    output: MULTIQC
    params: time="2:00:00", mem="1000m",
        outDir = MULTIQC_DIR,
        searchDir = DATA_DIR
    threads: 1
    conda: CONFIG_DIR + "multiqc_env1.yaml"
    shell:
        """
        multiqc -o {params.outDir} \
        {params.searchDir}
        """
