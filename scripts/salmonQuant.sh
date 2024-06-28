#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

# load salmon in /mnt/storage/apps
export PATH=/mnt/storage/apps/salmon-1.10.1/bin:$PATH


reads=(${snakemake_input[reads]})
r1="${reads[0]}"
r2="${reads[1]}"

GENOMEDIR="${snakemake_input[index]}"
OUTDIR=$(dirname "${snakemake_output[0]}")

# Create output directory if it does not exist
mkdir -p $OUTDIR

salmon quant -i $GENOMEDIR \
-l A \
-p 30 \
-1 "${r1}" \
-2 "${r2}" \
-o "${OUTDIR}/" \
--validateMappings \
--gcBias \
--seqBias
