#FROM gencommand_image:v1.0.0
FROM ubuntu:22.04

# Set up working directory
WORKDIR /pipeline

# Install basic dependencies
RUN apt-get update && apt-get install -y wget bzip2 && apt-get clean

# Install Miniconda (Linux x86_64)
ENV CONDA_DIR=/opt/conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p $CONDA_DIR && \
    rm miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# Install mamba and bio tools
RUN conda install -y -n base -c conda-forge mamba && \
    mamba create -y -n rnaseq-env -c conda-forge -c bioconda -c defaults \
        fastqc \
        multiqc \
        fastp \
        star \
        subread \
        r-base \
        r-ggplot2 \
        bioconductor-deseq2 \
        r-pheatmap

# Activate environment by default
SHELL ["conda", "run", "-n", "rnaseq-env", "/bin/bash", "-c"]

CMD ["/bin/bash"]
