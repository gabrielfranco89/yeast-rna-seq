#!/bin/bash
set -e

echo "📥 Downloading yeast reference genome and GTF..."

mkdir -p data/genome
cd data/genome

# Download the genome FASTA
wget https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

# Download the GTF annotation file
wget https://ftp.ensembl.org/pub/release-109/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz

# Unzip both files
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz

echo "✅ Files downloaded and unzipped."

# Create STAR index
mkdir -p star_index

echo "🧠 Building STAR genome index..."
STAR --runMode genomeGenerate \
  --genomeDir star_index \
  --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.109.gtf \
  --sjdbOverhang 75 \
  --runThreadN 4

echo "✅ STAR genome index built successfully."
cd -
