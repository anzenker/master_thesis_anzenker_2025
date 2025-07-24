#  Transcriptome Assembly and Annotation with ONT dRNA-seq data

###### ================================================================================
###### OVERVIEW
+ 1. [ This Repository ](#rep)
+ 2. [ Hardware and Operating System ](#OS)
+ 4. [ Guide to assemble and annotate a Transcriptome from ONT dRNA-seq data ](#guide)
    * 4.1 [ Basecalling & Preprocessing ](#prepros)
    * 4.2 [ nanoTome Pipeline ](#nanotome)
    * 4.3 [ Additional Annoation ](#addanno)
+ 5. [ Some file formats explained ](#file-formats)
+ 6. [ Some extra commands ](#extra)
+ 7. [ Abbreviations ](#abbrev)
###### ================================================================================

<a name="rep"></a>
## 1. This Repository
This repository offers a guide to assemble a Transcriptome from Oxford Nanopore Technologies (ONT) direct RNA (dRNA) sequencing data with a reference genome as a guide.

This guide is separated into 3 Parts:
1. Basecalling & Preprocessing of the Raw Sequencing Data
2. nanoTome pipeline: Transcriptome Assembly & Annotation
3. Additional Annotation

<a name="OS"></a>
## 2. Hardware and Operating System 
All steps except dorado basecalling were successfully executed on Linux-based systems with a x86_64 architecture.

Dorado basecalling was executed on Linux-based systems (x86_64) equipped either with NVIDIA RTX 6000 Ada Generation or NVIDIA Tesla V100-DGXS-32GB GPUs.

<a name="install"></a>
## 3. Installation of required software 

### Software which must be installed manually:
- [4.1 General Software Requirements](#gensof)
- [4.2 Preprocessing](#prepro)
- UniProt Annotation via blastp

For the rest there are two options of installing the required software:

### Option 1: Manual Installation
You can install each tool individually.
For each tool, a link to the official installation instructions is provided.
- [4.3 Transcriptome Reconstruction](#trans-recon)
- [4.4 Transcriptome Annotation](#trans-anno)

For this option go to [README_manual_software_installation.md](README_manual_software_installation.md)

### Option 2: Use the available Docker Image (recommended for main analysis):
For steps [4.3 Transcriptome Reconstruction](#trans-recon) and [4.4 Transcriptome Annotation](#trans-anno), a prebuilt Docker Image is available. It contains all necessary tools (except UniProt annotation).
This Docker image can be used together with the provided Nextflow pipeline for easy and reproducible analysis.
A description of how this is done is provided in [`README_Master_Thesis_commands_used.md`](README_Master_Thesis_commands_used.md).

<a name="guide"></a>
## 4. Guide to assemble and annotate a Transcriptome from ONT dRNA-seq data

<a name="prepros"></a>
## 4.1. Basecalling & Preprocessing
xxx

<a name="nanotome"></a>
## 4.2. nanoTome Pipeline

The nanoTome pipeline is designed to assemble transcriptomes from Oxford Nanopore Technologies (ONT) direct RNA sequencing (dRNA-seq) data, using a reference genome as a guide.

It combines steps for reconstructing the transcriptome, evaluating its completeness with BUSCO Vertebrata, and performing functional annotation with eggNOG. By this, the workflow provides a reproducible and user-friendly solution that saves time for assessing the quality and completeness of dRNA-seq data for further downstream analyses.

The nanoTome pipeline originates from a Master thesis project studying hybridization effects in the parthenogenetic species Aspidoscelis neomexicanus and its sexual parental species, A. marmoratus and A. arizonae. Since assembling a single transcriptome involves more than one time-consuming step, the aim was to combine and automate most of the process for efficient analysis of multiple tissue samples and species. The pipeline enables quicker and standardized transcriptome assembly, preparing the data for subsequent analyses.

<a name="pipe_install"></a>
#### 4.2.1 Installation Requirements
These tools are needed to set up the environment and run the pipeline.
- **[miniforge3](https://github.com/conda-forge/miniforge)**: Used to create environments with conda and mamba for easier package managment and software installation.
- **[nextflow](https://www.nextflow.io/docs/latest/install.html)**: Workflow manager to run the analysis pipeline in a reproducible way.
- **[docker](https://docs.docker.com/engine/install/ubuntu/)**: Allows to run software in containers. Required for using the provided Docker Image.
- **Java [conda_java20_env.yml](conda_java20_env.yml)**: OpenJDK version 20. Required to run Nextflow pipeline. Installed in a conda environment.


**Download and prepare necessary databases for the pipeline:**

- **EggNOG** (~13G & ~9G & ~7G)
```
eggnog link to database:
http://eggnog5.embl.de/download/emapperdb-5.0.2/

#download the databases into the 'bin' folder
cd bin/

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
gunzip eggnog.db

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
gunzip eggnog_proteins.dmnd

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
tar –xvzf eggnog.taxa.tar.gz
```
- **BUSCO** (~ 529M)
```
busco link to vertebrata odb10 database
https://busco-data.ezlab.org/v5/data/lineages/vertebrata_odb10.2024-01-08.tar.gz

#download the databases into the 'bin' folder
cd bin/

wget https://busco-data.ezlab.org/v5/data/lineages/vertebrata_odb10.2024-01-08.tar.gz
tar -xzf vertebrata_odb10.2024-01-08.tar.gz
```

<a name="pipe_run"></a>
#### 4.2.2 Run the nanoTome pipeline
The nextflow pipeline runs on a Docker image by default.

```
# help message
nextflow run main.nf -help

#run
nextflow run main.nf --raw_reads raw_reads.fastq --genome genome.fa --threads NO_THREADS
```

*Optional parameters 
- --outdir NAME     (defines the anme of the output directory) 
- --color HEX_CODE  (defines the plot color for the output plots, requires a hex color code e.g. #688e26)
- skip_eggnog true  (runs the workflow without eggNOG annotation)
- skip_busco true   (runs the workflow without busco analysis)
- skip_orf          (runs the workflow without orf prediction and without eggNOG annotation)
- no_plots          (no outout plots are generated from the workflow)
*

<a name="pipe_dir"></a>
#### 4.2.3 nanoTome results directory
```
.
└──  RESULTS
│   ├── 1_minimap2_output
│       └── ...
│   ├── 2_stringtie2_transcriptome
│       └── ...
│   ├── 3_gffread_transcriptome
│       └── ...
│   ├── 4_canonical_transcriptome
│       └── ...
│   ├── 5_frame_selection
│       └── ...
│   ├── 6_busco_vertebrata_completeness
│       ├── busco_output_xxx_transcripts
│           └── ...
│       ├── busco_output_xxx_transcripts_caninical
│           └── ...
│   ├── 7_eggnog_annotation
│           └── ...
│   ├── 8_plots
│           └── ...
```
<a name="manual"></a>
#### 4.2.4 Code implemented into the nanoTome pipeline
All code implemented intot the pipeline can be found in [README_commands_implemented_in_pipeline.md](README_commands_implemented_in_pipeline.md) for manual execution.

<a name="addanno"></a>
## 4.3. Additional Annotation
<a name="up"></a>
### 4.3.1 UniProt annotation 

### Download and make database locally
```
UniProt link to database:
https://www.uniprot.org/help/downloads

#download the databases into the 'bin' folder
cd bin/

# download database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzip downloaded fasta file
gzip -d uniprot_sprot.fasta.gz

# create blast database from fasta file (~310M)
makeblastdb -in bin/uniprot_sprot.fasta -dbtype prot

```

### Run the annotation
```
# blast search of transcripts against UniProt database
# choose max_target_seq 1 to only receive the top hit

blastp -query transdecoder_dir/longest_orfs.pep  \
    -db uniprot_sprot.fasta  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
```

<a name="file-formats"></a>
## 5. Some file formats explained
Some file formats used to analyse the data are explained in  [`README_file_formats.md`](README_file_formats.md).


<a name="extra"></a>
## 6. Some extra commands
Command to check the amount of bases in the .fasta or .fastq files:
```
seqkit stats *.fasta
```
Command to get the sequence ID out of .fasta `"^>"` file. This was useful to find out if the target species had en extra plasmid next to its chromosome or not.
```
grep "^>" dateiname.fasta
```
Command to get amount of reads from .fastq.gz `"^@"` file
```
zcat your_file.fastq.gz | grep -c "^@"
```

<a name="abbrev"></a>
## 7. Abbreviations
- ONT   Oxford Nanopore
- QC    Quality Control
- mRNA  messenger RNA
