#!/bin/bash
mkdir -p results/alignments
GENOME_DIR=data/genome/star_index

for fq1 in data/clean_fastq/*_1.trimmed.fastq.gz; do
  fq2="${fq1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
  base=$(basename "$fq1" _1.trimmed.fastq.gz)

  STAR \
    --genomeDir $GENOME_DIR \
    --readFilesIn "$fq1" "$fq2" \
    --readFilesCommand zcat \
    --outFileNamePrefix results/alignments/${base}_ \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 4
done

