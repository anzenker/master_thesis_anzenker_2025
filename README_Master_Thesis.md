#  Master Thesis - Transcriptome Analyses in Aspidoscelis neomexicanus

###### ================================================================================
###### OVERVIEW
+ 1. [ This Repository ](#rep)
+ 2. [ Hardware and Operating System ](#OS)
+ 3. [ Installation of required software  ](#install)
    * 3.1 [ General software requirements ](#gensof)
    * 3.2 [ Preprocessing ](#prepro)
    * 3.3 [ Transcriptome Reconstruction ](#trans-recon)
    * 3.4 [ Transcriptome Annotation ](#trans-anno)
+ 4. [ Guide to Assemble and Annotate a Transcriptome from ONT dRNA-seq data ](#guide)
+ 1. Preprocessing
+ 2. [ nanoTome Pipeline ](#nanotome)
+ 3. 

+ 5. [ Steps to recreate the results ](#steps)
+ 6. [ Some File Formats explained ](#file-formats)

+ 9. [ Some extra commands ](#extra)
+ 10. [ Abbreviations ](#abbrev)
###### ================================================================================

<a name="rep"></a>
## 1. This Repository
This repository offers a guide to assemble a Transcriptome from Oxfort Nanopore Technologies (ONT) direct RNA (dRNA) sequencing data with a reference genome as a guide.

This guide is separated into 3 Parts:
1. Basecalling & Preprocessing of the Raw Seqeuncing Data
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






<a name="nanotome"></a>
## 2. nanoTome Pipeline

The nanoTome pipeline is designed to assemble transcriptomes from Oxford Nanopore Technologies (ONT) direct RNA sequencing (dRNA-seq) data, using a reference genome as a guide.

It combines steps for reconstructing the transcriptome, evaluating its completeness with BUSCO Vertebrata, and performing functional annotation with eggNOG. By this, the workflow provides a reproducible and user-friendly solution that saves time for assessing the quality and completeness of dRNA-seq data for further downstream analyses.

The nanoTome pipeline originates from a Master thesis project studying hybridization effects in the parthenogenetic species Aspidoscelis neomexicanus and its sexual parental species, A. marmoratus and A. arizonae. Since assembling a single transcriptome involves more than one time-consuming step, the aim was to combine and automate most of the process for efficient analysis of multiple tissue samples and species. The pipeline enables quicker and standardized transcriptome assembly, preparing the data for subsequent analyses.

<a name="pipe_install"></a>
#### 2.1 Installation Requirements
These tools are needed to set up the environment and run the pipeline.
- **[miniforge3]()**: Used to create environments with conda and mamba for easier package managment and software installation.
- **[nextflow]()**: Workflow manager to run the analysis pipeline in a reproducible way.
- **[docker]()**: Allows to run software in containers. Required for using the provided Docker Image.
- **Java [conda_java20_env.yml](conda_java20_env.yml)**: OpenJDK version 20. Required to run Nextflow pipeline. Installed in a conda evironment.

<a name="pipe_run"></a>
#### 2.2 Run the nanoTome pipeline
The nextflow pipeline runs on a Docker image by default.

```
nextflow run main.nf --raw_reads raw_reads.fastq --genome genome.fa --threads NO_THREADS --outdir OPTIONAL --color OPTIONAL
```

*Optionally the paramters --outdir and --color can be given to the command. Be aware that --color requires a hex color code e.g. #688e26*

<a name="pipe_dir"></a>
#### 2.3 nanoTome results directory
```
.
└──  RESULTS
│   ├── **1_minimap2_output**
│       ├── xxx_mapped.bam
│       └── xxx_mapped.bam.bai
│   ├── 2_stringtie2_transcriptome
│       └── stringtie2_xxx_transcripts.gtf
│   ├── 3_gffread_transcriptome
│       └── stringtie2_xxx_transcripts.fasta
│   ├── 5_canonical_transcriptome
│       ├── stringtie2_xxx_transcripts_canonical.fasta
│       ├── stringtie2_xxx_transcripts_canonical_ids.txt
│       └── transcript_ids_and_coverage.tsv
│   ├── 6_frame_selection
│       └── stringtie2_xxx.transdecoder_dir
│           └── longest_orfs.pep
│   ├── 7_busco_vertebrata_completeness
│   ├── 8_plots
│       ├── busco_output_xxx_transcripts
│       ├── logs/
│           └── ...
│       ├── run_vertebrata_odb10/
│           ├── full_table.tsv
│           ├── short_summary.txt
│           └── ...
│       └── tmp/
│           └── ...
│       ├── busco_output_xxx_transcripts_canonical
│           ├── ... (same as other busco folder)
│   ├── 8_plots
│       ├── no_all_transcripts_vs_canonical_transcripts.png
│       ├── stringtie2_AC73_Amarm_Q10_noR_mapped_transcripts_canonical_length_distribution.png
│       ├── stringtie2_xxx_isoform_per_gene_barplot.png
│       ├── stringtie2_xxx_isoform_per_gene_barplot.txt
│       └── stringtie2_AC73_Amarm_Q10_noR_mapped_transcripts_length_distribution.png
```
...



## 4.5 Comparative Analyses
????

<a name="steps"></a>
## 5. Steps to recreate the results
The steps are explained in [`README_Master_Thesis_commands_used.md`](README_Master_Thesis_commands_used.md).

<a name="file-formats"></a>
## 6. Some file formats explained
Some file formats used to analyse the data are explained in  [`README_file_formats.md`](README_file_formats.md).



<a name="extra"></a>
## 9. Some extra commands
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
Command to split .fastq.gz into specified no. of parts (-p)
```
seqkit split2 reads.fq.gz-p 2 -O out -f --by-part-prefix "x_r{read}_"
    
```

<a name="abbrev"></a>
## 9. Abbreviations
- ONT   Oxford Nanopore
- QC    Quality Control
- mRNA  messenger RNA
