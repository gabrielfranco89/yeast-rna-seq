#!/bin/bash
mkdir -p results/counts

GTF=data/genome/Saccharomyces_cerevisiae.R64-1-1.109.gtf

featureCounts \
  -a "$GTF" \
  -o results/counts/count_matrix.txt \
  -T 4 \
  -p \
  -g gene_id \
  -t exon \
  -s 0 \
  results/alignments/*Aligned.sortedByCoord.out.bam
