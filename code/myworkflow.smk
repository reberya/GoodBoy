################################################################################
# TITLE:                    2020 RNAseq pipeline
# AUTHOR:                   Ryan Rebernick
# DATE LAST MODIFIED:       03/30/2020
#
# FUNCTION:                 FastQC
#
# USES:
#                           - python3.7-anaconda/2019.07
#                           - smk_env1
#                           - fastQC
#
# RUN: snakemake -s code/myworkflow.smk --cores 1
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
## samples, = glob_wildcards("data/sample_{sample}.txt")

##############################################################
# List of directories needed and end point files for analysis
##############################################################

FQC = expand("/home/rebernrj/testenv/output/fastqc/{sample}.R{read}_fastqc.html", sample=SAMPLES, read=READS)


#OLD: 127811_CGAATACG-TTACCGAC_S1_R1_001.fastq.gz
#NEW: 127811_CGAATACG-TTACCGAC_S1.R1.fastq.gz

##############################################################
# Rules
##############################################################


# end files required
rule all:
    input: FQC
    params: time="10:00:00"


# FastQC
rule fastqc:
    input:
        files = np.unique(expand("/home/rebernrj/testenv/data/FQ/{sample}.R{read}.fastq.gz", sample=SAMPLES, read=READS))
    output:
        reads = "/home/rebernrj/testenv/output/fastqc/{sample}.R{read}_fastqc.html",
    params: time="02:00:00", mem="4000m"
    threads: 2
    shell: """ \
    mkdir -p /home/rebernrj/testenv/output/fastqc; \
    fastqc -o /home/rebernrj/testenv/output/fastqc/ {input.files} \
    """
