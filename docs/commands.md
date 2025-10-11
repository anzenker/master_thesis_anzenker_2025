#  Commands Used

###### ================================================================================
###### OVERVIEW
* **1. [Preprocessing](#prepro)**
    * 1.1. [ dorado basecalling ](#dorado)
    * 1.2. [ ubam to fastq ](#file_format)
    * 1.3. [ Filtering RNA CS Control Sequences ](#rna_cs)
    * 1.4. [ pycoQC ](#pycoQC)
    * 1.5. [ Filtering based on PHRED-score ](#phred_filter)
    * 1.6. [ NanoComp ](#NanoComp)
    * 1.7. [ poly(A) tail trimming](#polya)
    * 1.8. [ Repetitive Pattern Detection ](#rep_pat)

* **2. [Nextflow Pipeline - Implemented Tool Commands](#ms_pipe)**
    * 2.1. [ minimap2 & samtools - Mapping the raw reads to the reference genome ](#minimap2)
    * 2.2. [ stringtie2 - Reference-guided Transcriptome Assembly ](#stringtie2)
    * 2.3. [ gffread - Reading the transcript sequences from the Reference Genome  ](#gffread)
    * 2.4. [ Choosing the Canonical Transcriptome ](#cano)
    * 2.5. [ BUSCO - Vertebrata Ortholog Assessment ](#busco)
    * 2.5. [ TransDecoder - ORF Prediction ](#transdec)
    * 2.6. [ eggNOG Annotation] (#eggnog)
    * 2.7. [ Custom Python Scripts - Visualization of some initial Results ](#python)
  
###### ================================================================================

<a name="prepro"></a>
## 1. Preprocessing
<a name="dorado"></a>
### 1.1. [dorado basecalling](https://github.com/nanoporetech/dorado)

Converts ONT raw .pod5 files into unaligned .ubam reads. Here, we will generate an unaligned BAM file (.ubam) as we are not directly mapping the sequencing data to a reference genome file.

- input: *.pod5
- output: .ubam

```
dorado basecaller --emit-moves --estimate-poly-a sup,m5C,inosine_m6A,pseU /pathTo/pod5 > xxx_basecalls.ubam
```

<a name="file_format"></a>
### 1.2. [samtools](https://www.htslib.org/doc/samtools.html) - ubam to fastq

Convert .ubam basecalled reads into .fastq format.

- input: .ubam
- output: .fastq

```
samtools fastq -t -@ 8 xxx_basecalls.ubam > xxx_basecalls.fastq
```

<a name="rna_cs"></a>
### 1.3. Filtering RNA CS Control Sequences

Remove ONT spike-in RNA CS control reads from the dataset using minimap2 & seqkit.

- input: raw_reads.fasta
- output: read_ids.txt, filtered.fastq/.fasta

**(1) Download the Sequences** 
Download the RNA CS sequences. Follow the link: [RNA CS (RNC)](https://nanoporetech.com/support/library-prep/RNA-and-cDNA/what-is-rna-cs-rcs)

**(2) Map the CS sequences against your raw sequencing data [minimap2](https://github.com/lh3/minimap2)**
```
minimap2 -ax map-ont rna_cs_control.fasta your_raw_reads.fasta > cs_mapped.sam
```

**(3) Extract read IDs**
The first column in a SAM file contains the read ID, and the fifth column contains the mapping quality (MAPQ). This command checks the first column to ensure it is not a header row (starts with @) and checks each 5th row to see if it is greater than 0. From all rows that comply with these conditions, the read ID is written to a TXT file.

```
awk '$1 !~ /^@/ && $5 > 0 {print $1}' cs_mapped.sam > read_ids.txt
```

**(4) Exclude read IDs from the raw sequencing data with [seqkit](https://github.com/shenwei356/seqkit)**
Use seqkit to keep all raw reads whose read ID matches any in the read_id.txt file.

```
seqkit grep -v -f read_ids.txt raw_reads.fastq -o raw_reads_filtered.fastq
```

<a name="#pycoQC"></a>
### 1.4. [pycoQC](https://a-slide.github.io/pycoQC/installation/)

Generate interactive QC metrics and plots from ONT sequencing summary files.

The sequencing summary file, which is best suited for this Quality Control report, is available with the raw sequencing data. It is generated from the sequencing run. It holds information on IDs, pore health, start time, and other details of the sequencing run and data.
A report can also be generated from .fastq or .fasta files, but it will show less information on the sequencing run.

```
pycoQC --summary file sequencing_summary.txt -o pycoQC_output.html 
```

Another helpful tool to use for this case is [NanoPlot](https://github.com/wdecoster/NanoPlot).

<a name="phred_filter"></a>
### 1.5. [chopper](https://github.com/wdecoster/chopper) - Filtering based on PHRED-score

Filters sequencing reads based on mean PHRED quality score thresholds or other parameters.

- input: raw_reads.fastq/.fasta
- output: filtered.fastq/.fasta

```
gunzip –c xxx.fastq.gz | chopper --quality 10 --threads 18 --input | gzip xxx_filteredQ10.fastq.gz 

```

<a name="nanoComp"></a>
### 1.6. [NanoComp](https://github.com/wdecoster/nanocomp)

Compare multiple ONT sequencing runs or datasets based on QC metrics.

```
# nanoComp report (compare fastq, fasta, summary)

NanoComp --summary xxx_summary.txt xxx_summary.txt xxx_summary.txt --names xxx_1 xxx_2 xxx_3 -o qc/nanocomp_xxx 
```

<a name="polya"></a>
### 1.7. Poly(A) tail trimming **[/python/scripts/polyA_pattern_detection.py]**
ONT direct RNA sequencing data includes poly(A) tails. Trimming them is not necessary when the raw reads are mapped to a reference genome, as they will be soft-clipped by minimap2. 
However, it was observed that trimming the poly(A) tails was helpful for the reconstruction of a de novo transcriptome without a reference genome (e.g., using RNAbloom [https://github.com/bcgsc/RNA-Bloom]), as differing tail lengths can influence the overlap between reads.
There are different available tools, such as dorado --estimate_polyA and Nanopolish, which are developed to detect poly(A) tail length, but poly(A) trimming was not found to be a standard usage tool. 
In this repository, I provide a custom script developed during the project to detect poly(A) tails by the AAA pattern and to trim the sequences based on this knowledge. It is important to note that this is not a validated tool! The Python script relies on the concept of detecting a AAAA pattern of a specified length (default: 10). Patterns are detected at the beginning, the middle, and the end of reads. Redundant detections of middle pattern detections, including those with start and end detection, are removed. Start and end patterns are simply trimmed off. Middle patterns are checked if a nucleotide sequence of < 20nt follows them, and then they are always considered as poly(A) end patterns. (This is because it was noticed that more than a few reads showed poly(A) middle patterns, which seemed to be tails, as a small number of non-A nucleotides followed them). I want to point out that this is not a validated script; however, it was tested using the best available knowledge with the data from this project, and in these cases, it produced reliable detection and trimming.

- input: raw_reads.fastq
- output: overview.tsv, raw_reads_noA.fasta
- 
```
python detect_and_trim_polyA.py raw_read.fastq
```

<a name="rep_pat"></a>
### 1.8. Repetitive Pattern Detection **[/python/scripts/repetitive_pattern_detection.py]**
ONT direct RNA sequencing generated reads from our sequencing run, showing repetitive patterns with lengths up to 500 kb. As these repetitive patterns do not contain valuable information for us, and it was observed that they do influence the mapping of the raw reads to the reference genome, it was decided to exclude reads that show a repetitive pattern. 
This is a simple script for detecting different repetitive patterns (e.g., 'AGAGA', 'CTCTC', etc.) of a minimum length (default=15) and extracting their read IDs into a .txt file. 
It should be noted that for our sequencing runs, excluding these reads showing a repetitive pattern was feasible, as only ~20% of the total of ~24 million reads were removed. This should be decided individually for each sequencing run.

- input: raw_reads.fastq / raw_reads.fasta
- output: read_ids.txt

```
# this can also be done on a FASTA file format
python repetitive_pattern_detection.py raw_reads.fastq

seqkit grep -v -f read_ids.txt raw_reads.fastq -o raw_reads_noR.fastq
```

<a name="ms_pipe"></a>

<a name="nf_pipe"></a>
## 2. Nextflow Pipeline - Implemented Tool Commands 
## Transcriptome Reconstruction & Annotation

As the reconstruction and analysis of biological data can be very individual, and also to provide the user with insights into how each step was used and in what order, all commands implemented into the NextFlow pipeline are listed in this README. With this, it is also easy to follow each step and/or use one or more steps individually if necessary. 

#### **Commands implemented into the pipeline:**

<a name="minimap2"></a>
## 1.1 [minimap2](https://lh3.github.io/minimap2/minimap2.html) & [samtools](https://www.htslib.org/doc/samtools.html) - Mapping the raw reads to the reference genome

```
# 1. minimap2 - mapping of the raw reads to the reference genome
minimap2 -y --MD -ax splice -uf -k14 -t 24 xxx_genome.fa raw_reads.fastq.gz | samtools sort -@ 24 -o xxx_genome_mapped.bam

samtools index xxx_genome_mapped.bam
```
<a name="stringtie22"></a>
## 1.2  [stringtie2](https://github.com/skovaka/stringtie2) - Reference-guided Transcriptome Assembly

```
# 2. stringtie2 - transcriptome reconstruction
stringtie -o stringtie2_xxx_transcripts.gtf -l NAME_PREFIX -L -p 24 xxx_genome_mapped.bam
```

<a name="gffread"></a>
## 2.3  [gffread](https://github.com/gpertea/gffread) - Reading the transcript sequences from the Reference Genome 

```
# 3. gffread - extract transcript sequences
gffread -w xxx_transcripts.fa -g xxx_genome.fa stringtie2_xxx_transcripts.gtf
```

<a name="cano"></a>
## 2.4 Choosing the Canonical Transcriptome

One representative (canonical) transcript per gene ID is chosen. Here, for each gene ID, the isoform with the highest nucleotide coverage (i.e., the most bases supported by aligned reads) was selected as the canonical transcript. This is done with a custom script.

```
# 4. Choose the canonical transcripts
python 4_choose_canonical_transcripts.py input.tsv output.txt

#extract these canonical IDs into a filtered transcriptome FASTA
seqkit grep -f output.txt transcriptome.fasta -o canonical_transcripts.fasta
```
_input_tsv_: TSV file holding the information of gene ID, transcript ID, and the stringtie2 calculated coverage value.
_output_txt_: output TXT file name/path, holding the canonical transcript IDs

<a name="busco"></a>
## 2.5  [BUSCO](https://busco.ezlab.org/busco_userguide.html) - Vertebrata Ortholog Assessment

Assesses completeness of the transcriptome based on conserved orthologs from the Vertebrata dataset.

```
# 4. BUSCO - Transcriptome Completeness Assessment
# offline:
busco -i xxx_transcripts.fa -l vertebrata_odb10 --download_path PATH/TO/busco_downloads/ -o OUTPUT_FILDER/ -m transcriptome --offline -c NO_THREADS
# with data download
busco -i xxx_transcripts.fa -l vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome -c NO_THREADS
```

<a name="transdec"></a>
## 2.6  [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) - ORF Prediction

Predicts coding regions (open reading frames, ORFs) from transcript sequences.

```
# 5. TransDecoder - predict Open Reading Frames (ORFs)
TransDecoder.LongOrfs -t xxx_transcripts.fa 
```

<a name="eggnog"></a>
## 2.7  [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13) Annotation

Annotates predicted protein sequences with orthology-based functional information.

(1) Download the database
For database download instructions, refer to [README_commands_used.md](README_commands_used.md)

(2) run the eggNOG command:
```
emapper.py  -m diamond --itype proteins -i xxx_transdecoder.pep -o NAME_PREFIX --data_dir /folder/with/databases
```

<a name="python"></a>
## 2.8  Custom Python Scripts - Visualization of some initial Results
- [plot](/python_scripts/4_plot_total_vs_canonical_transcript_count.py): count Total Transcriptome vs count Canonical Transcriptome
- [plot](/python_scripts/2_plot_isoform_per_gene.py): Isoforms per Gene
- [plot](/python_scripts/5_plot_orf_statistics.py): ORF Category distribution
- [plot](/python_scripts/6_busco_completeness_stacked_barplot.py): BUSCO Completeness
- [plot](/python_scripts/plot_pipeline_quality_overview.py): Overview Transcriptome Functional Assessment (Count Canonical Transcriptome, %with ORF, %with eggNOG Annotation, BUSCO Completeness)

