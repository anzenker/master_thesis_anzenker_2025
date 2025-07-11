FROM continuumio/miniconda3

# Set working directory
WORKDIR /testing

# Prepare conda
RUN conda init bash
RUN conda config --set auto_activate_base false
RUN conda install -n base -c conda-forge mamba
#RUN mamba install -c conda-forge -c bioconda busco=5.8.2
#RUN mamba install bioconda::seqkit bioconda::transdecoder
# Install dependencies via conda
RUN mamba install -y -c conda-forge -c bioconda \
    seqkit \
    transdecoder \
    python=3.10 \
    pip \
    biopython \
    pandas \
    bbmap \
    blast \
    augustus \
    metaeuk \
    prodigal \
    hmmer \
    sepp \
    r-base \
    r-ggplot2 \
    eggnog-mapper \
    blast \
 && conda clean -a

# Copy environment definition into container
#COPY environment_test.yml /tmp/environment_test.yml

# Create the environment and clean cache
#RUN conda env create -f /tmp/environment_test.yml && conda clean -a

# Set PATH so the environment is used
#ENV PATH /opt/conda/envs/busco_env/bin:$PATH
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pandas matplotlib Bio

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

######################
#install dependencies#
######################
RUN apt-get update -y && apt-get install -y \
    apt-utils \
    bzip2 \
    gcc \
    make \
    ncurses-dev \
    wget \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    ca-certificates \
    curl \
    g++ \
    gawk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
#***************************************************************************************************

# Install Java
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y openjdk-17-jdk-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
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
#seqkit 2.10.0#
################
ENV SEQKIT_VERSION=2.10.0
ENV SEQKIT_INSTALL_DIR=/opt/seqkit_f
# download & extract & install
WORKDIR /opt
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    rm seqkit_linux_amd64.tar.gz && \
    cp seqkit seqkit_f
# Add to PATH
ENV PATH="$PATH:${SEQKIT_INSTALL_DIR}"
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


