# Yeast RNA-Seq exercise

Repository to exercise a RNA-Seq pipeline, from public available data to DESeq2 analysis. 
Data comes from [Chen, Boning et al. (2025)](https://www.jbc.org/article/S0021-9258(25)00285-6/fulltext) study: _Spt6-Spn1 interaction is required for RNA Polymerase II association and precise nucleosome positioning along transcribed genes_.

For this study we used only RNA-Seq available and strings

- spt6_F249K, highly disruptive mutation
- spn1_R263D, mild disruptive mutation
- WT, wild type

## Reproducibility

You can run the scripts in the folder `scripts`, starting by the `download_fastq.sh` and then running them in order of their numbers.
It should be fine to run it in a PC with 16gb ram and 60Gb of free space. The whole process can take a couple hours, so be patient. 

The pipeline was run in a Docker container, which can be shared over request.

The code assumes the folder already exists, so check if you have them sorted before running.
