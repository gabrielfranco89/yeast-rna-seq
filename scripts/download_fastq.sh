#!/bin/bash
set -e

# Create target directory for FASTQ files
# mkdir -p data/raw_fastq

# List of SRR IDs (from GSE284959)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28495
# Files are TOO BIG
# Select only for for the same gene to just illustrata the pipeline

SRR_LIST=(
SRR31781631
SRR31781632
SRR31781633
#  SRR31781635
  # SRR31781636
 SRR31781640
  SRR31781641
)

# Download and convert using SRA Toolkit
for SRR in "${SRR_LIST[@]}"; do
  echo "🔽 Downloading and converting $SRR"

  # Download SRA file
  prefetch $SRR

  # Convert to FASTQ
  fasterq-dump $SRR --outdir data/raw_fastq --skip-technical  --threads 4
  
  # Compress
  echo "Compressing $SRR"
  gzip "data/raw_fastq/${SRR}_1.fastq"
  gzip "data/raw_fastq/${SRR}_2.fastq"
  
  # ToDO: Remove SRA file (space is not infinite)  

  echo "✅ Done: $SRR"
  echo "============================"
done
