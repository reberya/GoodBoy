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
################################################################################

from os.path import join
import os
import glob


##############################################################
# Globals
##############################################################

# samples/read names
READS = ["1", "2"]
SAMPLES = [os.path.basename(fname).split('_')[0] for fname in glob.glob('FQ/*fastq.gz')]
## samples, = glob_wildcards("data/sample_{sample}.txt")

##############################################################
# List of directories needed and end point files for analysis
##############################################################

DIRS = ['fastqc/', 'log/']
#FQC = expand("fastqc/{sample}.R{read}_val_{read}_fastqc.html", sample=SAMPLES, read=READS)
FQC = "data/TESTING.txt"

#127811_CGAATACG-TTACCGAC_S1_R1_001.fastq.gz

##############################################################
# Rules
##############################################################

localrules: dirs

# end files required
rule all:
    input:
        DIRS + FQC
    params: time="10:00:00"


# directories
rule dirs:
    output: DIRS
    params: time = "01:00:00"
    shell:  "mkdir -p "+' '.join(DIRS)

# FastQC
rule fastqc:
    output: FQC
    params: time= "01:00:00"
    shell: """ \
    echo "cat" > data/TESTING.txt
    """
