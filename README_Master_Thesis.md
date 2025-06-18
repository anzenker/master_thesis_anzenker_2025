#  Master Thesis - Transcriptome Analyses in Aspidoscelis neomexicanus

###### ================================================================================
###### OVERVIEW
+ 1. [ Task ](#task)
+ 2. [ Used data sets ](#data)
+ 3. [ Hardware and Operating System ](#OS)
+ 4. [ Installation of required software ](#install)
    * 4.1 [ miniforge3 ](#miniforge3)
    * 4.2 [ SRA toolkit ](#sra)
    * 4.3 [ rrwick Perfect bacterial-genome-tutorial (easy) ](#easy)
    * 4.4 [ pmenzel score-asseblies ](#sc)
+ 5. [ Steps to recreate the results ](#dir)
    * 5.1 [ rrwick Perfect bacterial-genome-tutorial (easy) ](#dirrr) 
        - 5.1.a. [ Directory Structure ](#dirrra)
        - 5.1.b. [ Scripts and Commands ](#dirrrb)
    * 5.2 [ pmenzel score-asseblies ](#dirpm)
+ 6. [ Short description of the assembly steps and their resulting files ](#assembly)
    * 6.1 [Read alignment with minimap2 and Samtools](#assembly1)
    * 6.2 [Read QC: Filtlong and fastp](#assembly2)
    * 6.3 [Consensus long-read assembler: Trycycler](#assembly3)
    * 6.4 [Long-read polisher: Medaka](#assembly4)
    * 6.5 [Short-read polishing: Polypolish using the BWA read aligner](#assembly5)
    * 6.6 [Short-read polishing: POLCA](#assembly6)
+ 7. [ Results presented with the tools badange, assemblytics and the snakemake-workflow score-assemblies ](#results)
+ 8. [ References](#ref)
+ 9. [ Some extra commands ](#extra)
###### ================================================================================

<a name="task"></a>
## 1. Scientific Objective 
The objective of this study is to investigate isoform diversity and parental gene
    contribution in the parthenogenetic lizard A. neomexicanus. Due to its hybrid
    genomic background, this species exhibits transcriptomic complexity that can
    be examined through long-read RNA sequencing.
    Using genome-guided transcriptome assembly with the use of a genome as-
    sembled from both parental species, the study aims to reveal patterns of parental
    contribution to gene expression. Integrated with functional annotation, isofroms
    and transcripts derived from one or both parental genomes will be examined on
    their biological function.
    All steps will be implemented in a reproducible Nextflow pipeline to enable
    future analyses of additional tissue samples or samples from the parental species. 

<a name="data"></a>
## 2. Pipeline
### 2.1 Preprocessing
- basecalling
- Quality control 
- Filtering & Trimming

### 2.2 Nextflow Pipeline
- minimap2
- stringtie2
- gffread
- python script - canonical transcriptome
- python scripts for plotting
- TransDecoder
- eggNOG
- interProScan
- UniProt

### 2.3 Comparative Analyses
- parental transcriptome vs hybrid transcriptome (no. of isoforms, no. of genes)
- parental gene contribution in hybrid transcriptome

<a name="data"></a>
## 3. Used data sets 
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
There are two options of installing the required software for this workflow.
First each software can be installed indivdually, for this each link to the original installation instructions i sprovided.
Second, the provided Dockerfile can be installed which already has all the necessary software installed in it.

 [README_install.md](README_install.md).
## 4.1 General software requirements
<a name="miniforge3"></a>
### 4.1.1 [miniforge3][1] [[1]] - conda and mamba are used to create environments for easier software installation
### 4.1.2 [nextflow][3] [[3]] - xxx
<a name="docker"></a>
### 4.1.3 [docker][4] [[4]] - xxx
<a name="xxx"></a>
### 4.1.4 [Java](conda_java20_env.yml) - conda environment with openjdk version "21.0.1" 2023-10-17, required to run Nextflow

## 4.2 Preprocessing 
<a name="dorado"></a>
### 4.2.1 [dorado][5] [[5]] - to basecall the raw sequencing data (.pod5 to .ubam) 
<a name="pycoQC"></a>
### 4.2.2 [pycoQC][6] [[6]] - generates metrics and interactive plots for ONT seqeuncing data
<a name="NanoPlot"></a>
### 4.2.3 [NanoPlot][7] [[7]] - 
<a name="NanoComp"></a>
### 4.2.4 [NanoComp][8] [[8]] - compare multiple runs of long read sequencing data and alignments

## 4.3 Transcriptome Reconstruction 
<a name="minimap2"></a>
### 4.3.1 [minimap2][9] [[9]] - map the raw reads to the genome
<a name="samtools"></a>
### 4.3.2 [samtools][10] [[10]] - work with .bam/.sam files
<a name="stringtie2"></a>
### 4.3.3 [stringtie2][11] [[11]] - reconstruct the transcriptome from the minimap2 mapping
<a name="gffread"></a>
### 4.3.4 [gffread][12] [[12]] - read the transcript seqeunces into a fasta file




<a name="dir"></a>
## 5. Steps to recreate the results

<a name="dirrr"></a>
### 5.1 [rrwick Perfect bacterial-genome-tutorial (easy)][4] [[4]]


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

[13]: ...

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




