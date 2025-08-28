FROM continuumio/miniconda3

LABEL maintainer="Annika Zenker <anzenker@students.uni-mainz.de"
LABEL description="Docker image for the ms-pipeline: long-read transcriptome analysis with Nanopore dRNA-seq data."
LABEL version="1.0"

# Prepare conda
RUN conda init bash
RUN conda config --set auto_activate_base false

# Install mamba
RUN conda install -n base -c conda-forge mamba

# Install dependencies via conda
RUN mamba install -y -c conda-forge -c bioconda \
    transdecoder \
    python=3.10 \
    pip \
    bbmap \
    blast \
    seqkit \
    augustus \
    metaeuk \
    prodigal \
    hmmer \
    sepp \
    eggnog-mapper \
 && conda clean -a

# Install python packages
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pandas matplotlib Bio

# Install dependencies
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    make \
    wget \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    ca-certificates \
    curl \
    g++ && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
#***************************************************************************************************

################
#Java#
################
# Install Java
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y openjdk-17-jdk-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
#***************************************************************************************************

################
#BUSCO#
################
# Clone BUSCO repo and install
WORKDIR /opt
RUN git clone https://gitlab.com/ezlab/busco.git
WORKDIR /opt/busco
RUN python -m pip install .

# Make sure busco is on PATH
ENV PATH="$PATH:/opt/busco/bin"
# Test BUSCO version
RUN busco --version
#***************************************************************************************************

################
#Samtools 1.3.1#
################
ENV SAMTOOLS_VERSION=1.22
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools
# download & extract & install
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar --bzip2 -xf samtools-1.22.tar.bz2 && \
    rm samtools-1.22.tar.bz2 && \
    mv samtools-1.22 samtools && \
    cd samtools && make
# Add to PATH
ENV PATH="$PATH:${SAMTOOLS_INSTALL_DIR}"
#***************************************************************************************************

################
#minimap2 2.29#
################
ENV MINIMAP2_VERSION=3.0.0
ENV MINIMAP2_INSTALL_DIR=/opt/minimap2
# download & extract & install
WORKDIR /opt
RUN wget https://github.com/lh3/minimap2/releases/download/v2.29/minimap2-2.29.tar.bz2 && \
    tar --bzip2 -xf minimap2-2.29.tar.bz2 && \
    rm minimap2-2.29.tar.bz2 && \
    mv minimap2-2.29 minimap2 && \
    cd minimap2 && make
# Add to PATH
#RUN echo $PATH
ENV PATH="$PATH:${MINIMAP2_INSTALL_DIR}"
#***************************************************************************************************


################
#stringtie 3.0.0#
################
ENV STRINGTIE_VERSION=3.0.0
ENV STRINGTIE_INSTALL_DIR=/opt/stringtie
# download & extract & install
WORKDIR /opt
RUN wget https://github.com/gpertea/stringtie/releases/download/v3.0.0/stringtie-3.0.0.tar.gz && \
    tar -xzf stringtie-3.0.0.tar.gz && \
    rm stringtie-3.0.0.tar.gz && \
    mv stringtie-3.0.0 stringtie && \
    cd stringtie && make
# Add to PATH
ENV PATH="$PATH:${STRINGTIE_INSTALL_DIR}"
#RUN echo $PATH
#***************************************************************************************************

################
#gffread 0.12.7#
################
ENV GFFREAD_VERSION=0.12.7
ENV GFFREAD_INSTALL_DIR=/opt/gffread
# download & extract & install
WORKDIR /opt
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.tar.gz && \
    tar -xzf gffread-0.12.7.tar.gz && \
    rm gffread-0.12.7.tar.gz && \
    mv gffread-0.12.7 gffread && \
    cd gffread && make
# Add to PATH
ENV PATH="$PATH:${GFFREAD_INSTALL_DIR}"
#***************************************************************************************************

################
#eggnog-mapper 2.1.13#
################
ENV EGGNOG_VERSION=2.1.13
ENV EGGNOG_INSTALL_DIR=/opt/eggnog
# download & extract & install
WORKDIR /opt
RUN wget https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.1.13.tar.gz && \
    tar -xzf 2.1.13.tar.gz && \
    rm 2.1.13.tar.gz && \
    cd eggnog-mapper-2.1.13 && python setup.py install
# Add to PATH
ENV PATH="$PATH:${EGGNOG_INSTALL_DIR}"    
#***************************************************************************************************


