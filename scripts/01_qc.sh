#!/bin/bash
mkdir -p results/qc_reports

fastqc data/raw_fastq/*.fastq -o results/qc_reports
multiqc results/qc_reports -o results/qc_reports
