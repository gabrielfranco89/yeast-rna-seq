#!/bin/bash
mkdir -p data/clean_fastq results/fastp_reports

for fq1 in data/raw_fastq/*_1.fastq.gz; do
  fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
  base=$(basename "$fq1" _1.fastq.gz)

  fastp \
    -i "$fq1" -I "$fq2" \
    -o data/clean_fastq/${base}_1.trimmed.fastq.gz \
    -O data/clean_fastq/${base}_2.trimmed.fastq.gz \
    --html results/fastp_reports/${base}_fastp.html \
    --json results/fastp_reports/${base}_fastp.json
done
