#  Commands Used

###### ================================================================================
###### OVERVIEW
+ 1. [Preprocessing](#prepro)
    * 1.1. [ dorado basecalling ](#dorado)
    * 1.2. [ ubam to fastq ](#file_format)
    * 1.3. [ pycoQC ](#pycoQC)
    * 1.4. [ Filtering based on PHRED-score ](#phred_filter)
    * 1.5. [ NanoComp ](#NanoComp)
    * 1.6. [ poly(A) tail trimming](#polya)
    * 1.7. [ Repetitive Pattern Detection ](#rep_pat)
+ 2. [Transcriptome Reconstruction & Annotation](#nf_pipe)
+ 3. [Analyses](#analyses)    

   
###### ================================================================================

<a name="prepro"></a>
## 1. Preprocessing
<a name="dorado"></a>
### 1.1. dorado basecalling [https://github.com/nanoporetech/dorado]
We will generate an unaligned BAM file (.ubam) as we are not directly mapping the sequencing data to a refence genome file.
- input: *.pod5
- output: .ubam
```
dorado basecaller --emit-moves --estimate-poly-a sup,m5C,inosine_m6A,pseU /pathTo/pod5 > xxx_basecalls.ubam
```

<a name="file_format"></a>
### 1.2. samtools - ubam to fastq [https://www.htslib.org/doc/samtools.html]
- input: .ubam
- output: .fastq
```
samtools fastq -t -@ 8 xxx_basecalls.ubam > xxx_basecalls.fastq
```

<a name="file_format"></a>
### 1.3. pycoQC [https://a-slide.github.io/pycoQC/installation/]
The sequencing summary file to best use for this Quality Control report is availabel with the raq seqeuncing data. It is generated from the sequencing run. It hold information of ids, pore healt, start time, etc.
A report can also be generated from .fastq or .fasta files but it will show less information on the sequencing run.
```
pycoQC --summary file sequencing_summary.txt -o pycoQC_output.html 
```

<a name="phred_filter"></a>
### 1.4. chopper - Filtering based on PHRED-score [https://github.com/wdecoster/chopper]

```
gunzip –c xxx.fastq.gz | chopper --quality 10 --threads 18 --input | gzip xxx_filteredQ10.fastq.gz 

```

<a name="nanoComp"></a>
### 1.5. NanoComp [https://github.com/wdecoster/nanocomp]

```
# nanoComp report (compare fastq, fasta, summary)

NanoComp --summary xxx_summary.txt xxx_summary.txt xxx_summary.txt --names xxx_1 xxx_2 xxx_3 -o qc/nanocomp_xxx 
```

<a name="polya"></a>
### 1.6. Poly(A) tail trimming **[CUSTOM SCRIPT]**
ONT direct RNA sequencing data includes poly(A) tails. Trimming them is not necessary when the raw reads are mapped to a reference genome, as they will be soft-clipped by minimap2. 
But as reconstruction a transcriptome for non-model organism still is quite experimental trimming the poly(A) tails was helpful in reconstruction a de novo transcriptome without a reference genome (f.e. with rnabloom [https://github.com/bcgsc/RNA-Bloom]) as differing tail length can influence the overlap between reads.
There are different available tools as dorado --estimate_polyA and Nanopolish which are developed to detect poly(A) tail length but poly(A) trimming was not found to be a common usage tool. 
In this repository I provide a custom script developed during the project to detect poly(A) tail by AAA pattern and to trim the sequences with this knowledge. It is important to not that this is not a validated tool. The python script relies on the idea of detecting a AAAA pattern of a specified length (default=10). Patterns are detected at the beginning, the middle, and at the end of reads. Redundant detections of middle pattern detections with the start and end detection are removed. Start and end patterns are simply trimmed off. Middle patterns are checked if they are followed by a nucleotide sequence of < 20nt and then they are always considered as poly(A) end patterns. (This ic because it was noticed that more than a few reads showed poly(A) middle patterns which seemd to be tails as they were follwed by a small number of not A nucleotides)
It is important to notice that this is not a validated script, but it was tested by the best knowledge with the data in this project and in these cases produced reliable detection and trimming.

- input: raw_reads.fastq
- output: overview.tsv, raw_reads_noA.fasta
```
python detect_and_trim_polyA.py raw_read.fastq
```

<a name="rep_pat"></a>
### 1.8. Repetitive Pattern Detection **[CUSTOM SCRIPT]**
ONT direct RNA sequencing generated in our sequenicng run reads showing repetitive pattern with lengths up to 500k bp. As these repetitive pattern do not hold any valuable infrmation for us and it was noticed that they do influence the mapping of the raw reads to the reference genome, it was decided to exclude reads showing a repetitive pattern.
This is a simple script for detecting different repetitive pattern (f.e. 'AGAGA', 'CTCTC', etc.) of a min length (default=15) and extracting their read ids into a .txt file. 
It should be noted that for our sequencing runs an exclusion of these reads showing a repetitive pattern was doable as only ~20% of a total of ~24M reads were removed. This should be decided individually for each seqeuncing run.

- input: raw_reads.fastq / raw_reads.fasta
- output: read_ids.txt

```
# this can also be done on fasta file format
python detect_and_trim_polyA.py raw_reads.fastq

seqkit grep -f read_ids.txt raw_reads.fastq -o raw_reads_noR.fastq
```

<a name="nf_pipe"></a>
## 2. Nextflow Pipeline - Transcriptome Reconstruction & Annotation
This pipeline is implemented in Nextflow (DSL2) using Docker for reproducibility. It supports transcriptome reconstruction and annotation from ONT dRNA-seq reads using reference-guided assembly and an initial annotation with EggNOG, InterProScan, and UniProt
.
The following directory structure and files are required to run the pipeline:
```
.
└──  SOME_FOLDER_NAME
│   ├── Dockerfile
│   ├── pipeline-test.nf
│   |── nextflow.config
│   ├── bin/
│       ├── eggnog.db
│       ├── eggnog_proteins.dmnd
│       ├── eggnog.taxa.db
│       ├── xxx UNIPROT DB ???
│       └── ...
|   ├── script/
│       ├── xxx
│       ├── xxx
│       └── ...
│   ├── pipeline_output
│       ├── xxx
│       ├── xxx
│       └── ...
│   ├──
│   └── ...
```
The pipeline integrates different tools providing an easy processing of raw data to a reconstructed transcriptome and its initial annotation. It is supposed to save time and provide the user with initial results and output files for further downstream analyses.

- input: raw_reads.fastq / raw_reads.fa, reference_genome.fa
- output: 
    - transcriptome.fa
    - canonical_transcriptome.fa
    - canonical_transritpome_annotated.gtf
    - Plot: busco_completeness.png, read_length_dist.png (before/after canonical reduction), isoforms_per_gene.png (before/after canonical reduction), orf_categories.png

There are two possible ways to use the pipeline
- either use it with a docker container, built from the provided Dockerfile
- or with each required software tools installed individually. 

Here I will only explain the usage with the docker container. The installation instructions for each software are linked in the **README_install.md**. 

```
# 1. docker - build docker container [https://docs.docker.com/engine/install/ubuntu/]
https://docs.docker.com/engine/install/ubuntu/

# Troubleshooting: I noticed that I could not run the nextflow pipeline as I has a Java version > v21 on my Oprating System. To solve this I created a simple conda environment with Java openjdk=20.0.2. For this I also provide the [conda_java20_env.yml] file from the environment I created. 

# either use the .yml file 
conda env create --name ENV_NAME --file=conda_java20_env.yml

# or create the environment without the -yml file
mamba env create --name ENV_NAME -c conda-forge openjdk=20.0.2

# 2. go into the working folder
cd SOME_FOLDER_NAME

# . nextflow - run the pipeline [https://www.nextflow.io/docs/latest/install.html]
nextflow run pipeline_test.nf --raw_reads raw_reads.fa --genome_file genome.fa --threads NO_THREADS --with-docker -c nextflow.config 


```
Note: If by any means the pipeline did not run in its entirety and the run want to be resumed after ficing the issue **-resume** can be added to the ```'nextflow run ...'``` command to continue were it has stopped int he previous run. 

To see the individual tool commands implemented into the pipeline go to [README_commands_implemented_in_pipeline.md].

<a name="analyses"></a>
## 3. Analyses
```
some analyses ...
```
