#!/bin/bash
set -e

# Create target directory for FASTQ files
# mkdir -p data/raw_fastq

# List of SRR IDs (from GSE284959)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28495
# Files are TOO BIG
# Select only for for the same gene to just illustrata the pipeline

SRR_LIST=(
#  SRR31781635
  SRR31781636
#  SRR31781640
  SRR31781641
)

# Download and convert using SRA Toolkit
for SRR in "${SRR_LIST[@]}"; do
  echo "🔽 Downloading and converting $SRR"

  # Download SRA file
  prefetch $SRR

  # Convert to gzipped FASTQ
  fasterq-dump $SRR --outdir data/raw_fastq --skip-technical  --threads 4

  echo "✅ Done: $SRR"
done
