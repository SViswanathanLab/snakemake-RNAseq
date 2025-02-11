#!/usr/bin/env bash
source /etc/profile.d/modules.sh
module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/samtools-1.9-gcc-6.5.0-b2repp2

exec 2> "${snakemake_log[0]}"  


read="${snakemake_input[bam_files]}"

THREADS=${snakemake_params[threads]}
Sorted_bam="${snakemake_output[0]}"
r1="${snakemake_output[1]}"
r2="${snakemake_output[2]}"

OUTDIR=$(dirname "${snakemake_output[0]}")
mkdir -p $OUTDIR # Create directory if it doesn't exist

samtools sort -n $read -o $Sorted_bam
samtools fastq -@ $THREADS $Sorted_bam -1 $r1 -2 $r2 -0 /dev/null -s /dev/null -n
