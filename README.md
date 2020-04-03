# GoodBoy
Basic paired end RNA-seq Analysis Pipeline for use on the University of Michigan Greatlakes server. 

## Getting Started

You will need access to the University of Michigan [Greatlakes server](https://arc-ts.umich.edu/greatlakes/) as well as a slurm account. Instructions for obtaining access to greatlakes and slurm are availible on the Greatlakes homepage.

The [Slurm User Guide](https://arc-ts.umich.edu/greatlakes/slurm-user-guide/) will be useful for job submission.


### Prerequisites
Once access is obtained to the Greatlakes server you may begin by cloning Goodboy into your directory of choice.
```
git clone <LINK_TO_GOODBOY>
```

### Installing

1. Load python from the cluster's existing modules. Version 3.7 with anaconda is preferred. If you choose a version without anaconda, then you will have to install [Anaconda](https://docs.anaconda.com/anaconda/install/) locally. 
```
module load python3.7-anaconda/2019.07
```
2. Create two conda environments from the availible yaml files.
```
conda env create --f Goodboy/config/python_env/multiqc_env1.yaml
```

```
conda env create --f Goodboy/config/python_env/smk_env1.yaml
```

### File modification
Several files must be modified to suit your specific installation of Goodboy. Additionally a reference genome for alignment and feature counting is required.

#### Goodboy.smk
There is a heading entitled 'User Modifiables' which contains 2 directories. 

1. DATA_DIR should be set equal to the parent folder of FQ/<your_fastq_files>.
```
DATA_DIR = 'path/to/your/folder/'
```

2. GOODBOY_DIR should be set equal to the location of your clone of Goodboy on the Greatlakes server.
```
GOODBOY_DIR = 'path/to/Goodboy/
```
#### submit_slurm.sh
The submit slurm.sh file should be modified for your specific account and mail-user.
```
#SBATCH --account=<your_slurm_account>
...
...
#SBATCH --mail-user=user@usermail.com

```

#### cluster.json
This file is located in Goodboy/config. Modify the account_slurm and email portions.
```
"account_slurm": "your_slurm_account",
...
...
"email": "user@usermail.com",
```

#### Reference genome
Make folder for reference genome in DATA_DIR. 
```
cd <DATA_DIR>
mkdir -p genome/GRCm38_vM24
cd genome/GRCm38_vM24
```

Get required files
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz
```

### Data organization
In order for Goodboy to run without error on the availible fastq files, the folders containing your data must be structured appropriately. Appropriate structure consists of a single folder 'FQ' which contains your .fastq.gz files in the following format: <SAMPLE>.R<#>.fastq.gz where SAMPLE is your sample name and # is the read number [1,2]. Example file names would be:
  
sample1.R1.fastq.gz
sample1.R2.fastq.gz

or...

14524_M.R1.fastq.gz
14524_M.R2.fastq.gz



## Running Goodboy
Once you have created the proper conda environments and structured your data appropriately, you may commence running Goodboy.  First navigate to Goodboy's directory. 

Next, you must do is activate the snakemake conda environment using:
```
conda activate smk_env1
```

Next, perform a dry run of snakemake to ensure there are no bugs due to above naming/structuring.
```
snakemake -s code/goodboy.smk --dryrun
```

If the dryrun completes without error, commence with slurm submission
```
sbatch code/submit_slurm.sh
```

You can check the status of your job using
```
squeue -u <your_username>
```



### Output
Goodboy will create an output directory in a folder within your DATA_DIR (see goodboy.smk) containing the output for each tool stored in a seperate folder. the output/featureCounts/counts.txt file may be used for gene expression with programs like [DESEQ2](https://github.com/mikelove/DESeq2) or [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

I recommend first starting with the MultiQC output to check sample quality and success of the run. It is located in output/multiqc/multiqc_report.html



## Built With

* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - FASTQC analysis
* [STAR](https://github.com/alexdobin/STAR) - Alignment of Fastqc files
* [SUBREAD](http://subread.sourceforge.net/) - Feature Counting
* [MULTIQC](https://multiqc.info/) - Overall QC (compiles samples)



## Acknowledgments

Shout out to user [Kelly Sovacool](https://github.com/kelly-sovacool) for her help setting up snakemake on greatlakes server initially. The initial configuration files and base snakemake file are based off of her mwe_hpc_snakemake repository. 
