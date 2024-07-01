#!/bin/bash

source /mnt/storage/apps/Mambaforge-23.1.0-1/bin/activate snakemake
#$ -S /bin/bash
#$ -q all.q
#$ -pe pvm 32
#$ -l h_vmem=64G
#$ -o $HOME/snakemake-RNAseq/joblogs/
#$ -e $HOME/snakemake-RNAseq/joblogs/

# Navigate to the directory containing the Snakefile
cd $HOME/snakemake-RNAseq

# Run Snakemake
snakemake --unlock
snakemake --cores all --latency-wait 60
