#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

reads=(${snakemake_input[reads]})  # don't double-quote this - we want word splitting
r1="${reads[0]}"
r2="${reads[1]}"

export PATH=/mnt/storage/apps/STAR/2.7.10a/bin:$PATH

GENOMEDIR="${snakemake_input[index]}"
GTFFILE="${snakemake_input[gtf]}"

OUTDIR=$(dirname "${snakemake_output[0]}")
sample_name=$(basename "$OUTDIR")
THREADS=${snakemake_params[threads]}

# Create output directory if it does not exist
mkdir -p $OUTDIR

if [[ "$r1" == *".gz"* ]]; then
  STAR --runThreadN $THREADS \
  --genomeDir $GENOMEDIR \
  --readFilesIn "${r1}" "${r2}" \
  --outFileNamePrefix "${OUTDIR}/${sample_name}_" \
  --readFilesCommand zcat \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --sjdbGTFfile $GTFFILE \
  --outReadsUnmapped Fastx \
  --outMultimapperOrder Random \
  --outWigType wiggle
else
  STAR --runThreadN $THREADS \
  --genomeDir $GENOMEDIR \
  --readFilesIn "${r1}" "${r2}" \
  --outFileNamePrefix "${OUTDIR}/${sample_name}_" \
  --outSAMattributes All \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM GeneCounts \
  --sjdbGTFfile $GTFFILE \
  --outReadsUnmapped Fastx \
  --outMultimapperOrder Random \
  --outWigType wiggle
fi
