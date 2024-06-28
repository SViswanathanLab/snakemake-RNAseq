#!/bin/bash

source /mnt/storage/apps/Mambaforge-23.1.0-1/bin/activate snakemake
#$ -S /bin/bash
#$ -q all.q
#$ -pe pvm 30
#$ -l m_mem_free=16000M
#$ -o $HOME/snakemake-RNAseq/joblogs/
#$ -e $HOME/snakemake-RNAseq/joblogs/

# Navigate to the directory containing the Snakefile
cd $HOME/snakemake-RNAseq

# Run Snakemake
snakemake --unlock
snakemake --cores 30 --latency-wait 60
