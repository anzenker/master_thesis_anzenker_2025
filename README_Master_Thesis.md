#  Master Thesis - Transcriptome Analyses in Aspidoscelis neomexicanus

###### ================================================================================
###### OVERVIEW
+ 1. [ Scientific Objective ](#sci-obj)
+ 2. [ Used data sets ](#data)
+ 3. [ Hardware and Operating System ](#OS)
+ 4. [ Installation of required software  ](#install)
    * 4.1 [ General software requirements ](#gensof)
    * 4.2 [ Preprocessing ](#prepro)
    * 4.3 [ Transcriptome Reconstruction ](#trans-recon)
    * 4.4 [ Transcriptome Annotation ](#trans-anno)
+ 5. [ Steps to recreate the results ](#steps)
+ 6. [ ... ](#...)
+ 8. [ References](#ref)
+ 9. [ Some extra commands ](#extra)
+ 10. [ Abbreviations ](#abbrev)
###### ================================================================================

<a name="sci-obj"></a>
## 1. Scientific Objective 
This study aims to investigate isoform diversity and parental gene
    contribution in the parthenogenetic lizard *A. neomexicanus*. Due to its hybrid
    genomic background, this species exhibits transcriptomic complexity, which can
    be examined through long-read RNA sequencing.
    
Using genome-guided transcriptome assembly with the use of a genome assembled from both parental species (*A. marmoratus* and *A. arizonae*), the goal is to reconstruct and compare transcriptomes of the parental species and *A. neomexicanus*, based on gene usage and isoform diversity. This is done using ONT direct RNA sequencing.

Transcriptome completeness will be assessed using BUSCO, and initial fucntional annotation will be performed (eggNOG, InterProScan, UniProt).

By comparing the transcriptomes of the parental species and their hybrid offspring, the study will provide insights into how hybridization affects isoform diversity and gene expression. 

All steps from Transcriptome Reconstruction and fuctional annotation (except InterProScan) will be implemented in a reproducible Nextflow pipeline to enable
    future analyses of additional tissue samples or samples from the parental species. 

<a name="data"></a>
## 2. Used data sets 
#### a. Oxford Nanopore direct RNA seqeuncing data 
    - Aspidoscelis neomexicanus - liver mRNA
    - SQK RNA-004
    - PromethION
#### b. Oxford Nanopore direct RNA seqeuncing data 
    - Aspidoscelis marmoratus - liver mRNA
    - SQK RNA-004
    - PromethION 
#### c. Genome Assembly 
    - Aspidoscelis marmoratus
    - Aspidoscelis arizonae
    
<a name="OS"></a>
## 3. Hardware and Operating System 
- Linux compsysgen2 5.15.0-140-generic #150-Ubuntu SMP; x86_64 x86_64 x86_64 GNU/Linux; Vendor ID: GenuineIntel, Model name: Intel(R) Xeon(R) Gold 5120 CPU @ 2.20GHz, CPU family: 6
- Linux compsysgen3 5.15.0-140-generic #150-Ubuntu SMP; x86_64 x86_64 x86_64 GNU/Linux; Vendor ID: GenuineIntel, Model name: Intel(R) Xeon(R) Gold 6438Y+, CPU family: 6
- Linux Olymp 5.4.0-212-generic #232-Ubuntu SMP, x86_64 x86_64 x86_64 GNU/Linux, Vendor ID: GenuineIntel, Model name: Intel(R) Xeon(R) CPU E5-2698 v4 @ 2.20GHz, CPU family: 6

<a name="install"></a>
## 4. Installation of required software 

### Software which must be installed manually:
- [4.1 General Software Requirements](#gensof)
- [4.2 Preprocessing](#prepro)
- [4.3.x InterProScan](#interproscan)

For the rest there are two options of installing the required software:

### Option 1: Manual Installation
You can install each tool individually.
For each tool, a link to the official installation instructions is provided.
- [4.3 Transcriptome Reconstruction](#trans-recon)
- [4.4 Transcriptome Annotation](#trans-anno)

### Option 2: Docker Image (recommended for main analysis):
For steps [4.3 Transcriptome Reconstruction](#trans-recon) and [4.4 Transcriptome Annotation](#trans-anno), a prebuilt Dockerfile is available. It contains all necessary tools (except InterProScan).
This Docker image can be used together with the provided Nextflow pipeline for easy and reproducible analysis.
A description of how this is done is provided in [`README_Master_Thesis_commands_used.md`](README_Master_Thesis_commands_used.md).

 <a name="gensof"></a>
## 4.1 General software requirements
These tools are needed to set up the environment and run the pipeline.
<a name="miniforge3"></a>
### 4.1.1 [miniforge3][1] [[1]]
Used to create environments with conda and mamba for easier package managment and software installation.
### 4.1.2 [nextflow][3] [[3]]
Workflow manager to run the analysis pipeline in a reproducible way.
<a name="docker"></a>
### 4.1.3 [docker][4] [[4]]
Allows to run software in containers. Required for using the provided Docker Image.
<a name="xxx"></a>
### 4.1.4 [Java](conda_java20_env.yml)
OpenJDK version 20. Required to run Nextflow pipeline. Installed in a conda evironment.

<a name="prepro"></a>
## 4.2 Preprocessing 
These tools are used before transcriptome reconstruction, to basecall, evaluate and refine the raw sequencing data.
<a name="dorado"></a>
### 4.2.1 [dorado][5] [[5]]
Used to convert raw ONT .pod5 files into .ubam basecalled reads.
<a name="pycoQC"></a>
### 4.2.2 [pycoQC][6] [[6]]
Generates metrics and interactive plots for ONT seqeuncing data.
<a name="NanoPlot"></a>
### 4.2.3 [NanoPlot][7] [[7]]
Generates metrics and plots for ONT seqeuncing data.
<a name="NanoComp"></a>
### 4.2.4 [NanoComp][8] [[8]]
Compares multiple runs of long read sequencing data and alignments based on QC metrics and plots.

<a name="trans-recon"></a>
## 4.3 Transcriptome Reconstruction 
<a name="minimap2"></a>
### 4.3.1 [minimap2][9] [[9]] - map the raw reads to the genome
<a name="samtools"></a>
### 4.3.2 [samtools][10] [[10]] - work with .bam/.sam files
<a name="stringtie2"></a>
### 4.3.3 [stringtie2][11] [[11]] - reconstruct the transcriptome from the minimap2 mapping
<a name="gffread"></a>
### 4.3.4 [gffread][12] [[12]] - read the transcript seqeunces into a fasta file

<a name="trans-anno"></a>
## 4.4 Transcriptome Annotation
<a name="eggNOG"></a>
### 4.3.1 [eggNOG][13] [[13]] - xxx
<a name="interproscan"></a>
### 4.3.2 [InterProScan][14] [[14]] - xxx
<a name="uniprot"></a>
### 4.3.3 [UniProt][15] [[15]] - xxx

## 4.5 Comparative Analyses
????

<a name="steps"></a>
## 5. Steps to recreate the results
The steps are explained in [README_Master_Thesis_commands_used.d].

<a name="dirrra"></a>
#### 5.1.a. Directory Structure
Below is the directory structure which will be generated by the Nextflow pipeline. 

```
.
└──  xxx
│   ├── xxx
│       ├── xxx
│       ├── xxx
│       ├── xxxx
│       └── xxx
│           ├── ...
│           ├── xxx
│               ├── ...
│               ├── xxx
│               └── ...
│           └── ...
│   ├── xxx
│   ├── xxx
│   └── xxx
```
...


<a name="ref"></a>
## 9. References
[1]: https://github.com/conda-forge/miniforge
[2]: https://www.nextflow.io/docs/latest/install.html
[3]: https://docs.docker.com/engine/install/ubuntu/
[4]: conda_java20_env.yml

[5]: https://github.com/nanoporetech/dorado
[6]: https://a-slide.github.io/pycoQC/installation/
[7]: https://github.com/wdecoster/NanoPlot
[8]: https://github.com/wdecoster/nanocomp

[9]: https://lh3.github.io/minimap2/minimap2.html "minimap2 user guide"
[10]: https://www.htslib.org/doc/samtools.html "samtools user guide"
[11]: https://github.com/skovaka/stringtie2
[12]: https://github.com/gpertea/gffread

[13]: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13
[14]: https://interproscan-docs.readthedocs.io/en/v5/
[15]: UniProt

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
