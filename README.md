#  Transcriptome Assembly and Annotation with ONT dRNA-seq data

###### ================================================================================
###### OVERVIEW
## 1. [ This Repository ](#rep)
## 2. [ Hardware and Operating System ](#OS)
## 3. [ Guide to assemble and annotate a Transcriptome from ONT dRNA-seq data ](#guide)
###    + 3.1 [ Basecalling & Preprocessing ](#prepros)
###    + 3.2 [ ms-pipeline (Master Thesis Pipeline) (ms-pipeline) ](#ms-pipeline)
####      + 3.2.1 [ Installation Requirements (install) ](#install)
####      + 3.2.2 [ Run the ms-pipeline (run) ](#run)
####      + 3.2.3 [ ms-pipeline flowchart overview (flowchart) ](#flowchart)
####      + 3.2.4 [ ms-pipeline results directory (dict) ](#dict)
####      + 3.2.5 [ Code implemented into the nanoTome pipeline (code) ](#code)
###    + 3.3 [ Additional Annoation with BLAST (addanno) ](#addanno)
####      + [ 3.3.1 UniProt annotation (uniprot) ] (#uniprot)
## 4. [ Some file formats explained ](#file-formats)
## 5. [ Links to software tools for manual installation (soft)](#soft)
## 6. [ Some extra commands ](#extra)
## 7. [ Abbreviations ](#abbrev)
###### ================================================================================

<a name="rep"></a>
## 1. This Repository
This repository provides a step-by-step guide for assembling and annotating a transcriptome from Oxford Nanopore direct RNA sequencing (dRNA-seq) data using a reference genome. It includes an initial functional assessment of the transcriptome, laying the foundation for reproducible analyses across different tissues, species, and replicates.

This guide is separated into 3 Parts:
1. Basecalling & Preprocessing of the Raw Sequencing Data
2. ms-pipeline: Transcriptome Assembly & Functional Assessment
3. Additional Annotation via BLAST/UniProt

<a name="OS"></a>
## 2. Hardware and Operating System 
All steps except dorado basecalling were successfully executed on Linux-based systems with an x86_64 architecture.

Dorado basecalling was executed on Linux-based systems (x86_64) equipped either with NVIDIA RTX 6000 Ada Generation or NVIDIA Tesla V100-DGXS-32GB GPUs.

<a name="guide"></a>
## 3. Guide to assemble and annotate a Transcriptome from ONT dRNA-seq data

<a name="prepros"></a>
## 3.1. Basecalling & Preprocessing
Basecalling and Preprocessing can be accomplished in different approaches. A useful approach for this is described in detail in [commands.md](/docs/commands.md).

<a name="ms-pipeline"></a>
## 3.2. ms-pipeline (Master Thesis Pipeline)

The ms-pipeline is designed to assemble transcriptomes from Oxford Nanopore Technologies (ONT) direct RNA sequencing (dRNA-seq) data, using a reference genome as a guide.

It combines steps for reconstructing the transcriptome (stringtie2), evaluating its completeness with BUSCO Vertebrata, predicting Open Reading Frames (TransDecoder) and performing functional annotation with eggNOG. By this, the workflow provides a reproducible and user-friendly solution that saves time for assessing the quality and completeness of the assembled Transcriptome from dRNA-seq data for further downstream analyses.

The ms-pipeline originates from a Master Thesis project studying hybridization effects in the parthenogenetic species *Aspidoscelis (A.) neomexicanus* and its sexual parental species, *A. marmoratus* and *A. arizonae*. Since assembling a single transcriptome involves more than one time-consuming step, the aim was to combine and automate most of the process for efficient analysis of multiple tissue samples and species. The pipeline enables quicker and standardized transcriptome assembly and its functional assessment, preparing the data for subsequent analyses.

<a name="install"></a>
#### 3.2.1 Installation Requirements
These tools are needed to set up the environment and run the pipeline.

- **[nextflow](https://www.nextflow.io/docs/latest/install.html)**: Workflow manager to run the analysis pipeline in a reproducible way.
- **[docker](https://docs.docker.com/engine/install/ubuntu/)**: Allows to run software in containers. Required for using the provided Docker Image [anzenker/ms-pipeline](https://hub.docker.com/repository/docker/anzenker/ms-pipeline/).

**Download and prepare necessary databases for the pipeline:**
The databases only need to be downloaded if these steps are to be executed with the workflow. The pipeline resolves the local path `/bin` at runtime

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

<a name="run"></a>
#### 3.2.2 Run the ms-pipeline
The nextflow pipeline runs on a Docker image by default.

```
# help message
nextflow run anzenker/master_thesis_anzenker_2025 -r main --help

#run
nextflow run anzenker/master_thesis_anzenker_2025 -r main --raw_reads raw_reads.fastq --genome genome.fa --threads NO_THREADS
```
**Quickstart - Run pipeline without eggNOG & BUSCO database download**
Information: This quick run will use a small test dataset from another repository.
```
nextflow run anzenker/master_thesis_anzenker_2025/ -r main -entry test -profile test
```
(runtime ~2 min, runs without BUSCO/eggNOG)

*Optional parameters 
- --outdir NAME     (defines the name of the output directory) 
- --color HEX_CODE  (defines the plot color for the output plots, requires a hex color code e.g. #688e26)
- --skip_eggnog true  (runs the workflow without eggNOG annotation)
- --skip_busco true   (runs the workflow without busco analysis)
- --skip_orf          (runs the workflow without orf prediction and without eggNOG annotation)
- --no_plots          (no output plots are generated from the workflow)
*

<a name="flowchart"></a>
#### 3.2.3 ms-pipeline flowchart overview
![ms_pipeline_flowchart.png](/images/ms_pipeline_flowchart.png)

An exemplary description of the output files and output plots can be found here: [exemplary_results_output.md](/docs/exemplary_results_output.md)

<a name="dict"></a>
#### 3.2.4 ms-pipeline results directory
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
│       ├── busco_output_xxx_transcripts_canonical
│           └── ...
│   ├── 7_eggnog_annotation
│           └── ...
│   ├── 8_plots
│           └── ...
```

<a name="code"></a>
#### 3.2.5 Code implemented into the nanoTome pipeline
All code implemented into the pipeline can be found in [commands.md](/docs/commands.md) for manual execution.

<a name="addanno"></a>
## 3.3. Additional Annotation
To generate a broader functional annotation of the Transcriptome.
<a name="uniprot"></a>
### 3.3.1 UniProt annotation 

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
## 4. Some file formats explained
Some file formats used to analyse the data are explained in  [file_formats.md](/docs/file_formats.md).

<a name="soft"></a>
## 5. Links to software tools for manual installation
[software_links.md](/docs/software_links.md)

<a name="extra"></a>
## 6. Some extra commands
Command to check stats in a .fasta or .fastq file:
```
seqkit stats *.fasta
```

<a name="abbrev"></a>
## 7. Abbreviations
- ONT   Oxford Nanopore
- QC    Quality Control
- mRNA  messenger RNA
