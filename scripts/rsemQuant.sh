#!/usr/bin/env bash
source /etc/profile.d/modules.sh
module load /mnt/storage/spack/share/spack/modules/linux-ubuntu16.04-x86_64/rsem-1.3.1-gcc-5.4.0-ml3p6ok

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

count_path="${snakemake_input[counts]}"
ref_path="rsem_ref/human_gencode"
count_file=$(basename "$count_path")
sample_name="${count_file%_Aligned.toTranscriptome.out.bam}"

OUTDIR=$(dirname "${snakemake_output[0]}")

mkdir -p $OUTDIR # create output directory if it does not exist

rsem-calculate-expression -p 30 \
--bam \
--seed 12345 \
--paired-end \
--no-bam-output \
"${count_path}" "${ref_path}" "${OUTDIR}/${sample_name}"
