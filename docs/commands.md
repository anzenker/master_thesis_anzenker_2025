#  Commands Used

###### ================================================================================
###### OVERVIEW
+ 1. [Preprocessing](#prepro)
    * 1.1. [ dorado basecalling ](#dorado)
    * 1.2. [ ubam to fastq ](#file_format)
    * 1.3. [ Filtering RNA CS Control Sequences ](#rna_cs)
    * 1.4. [ pycoQC ](#pycoQC)
    * 1.5. [ Filtering based on PHRED-score ](#phred_filter)
    * 1.6. [ NanoComp ](#NanoComp)
    * 1.7. [ poly(A) tail trimming](#polya)
    * 1.8. [ Repetitive Pattern Detection ](#rep_pat)

+ 2. [Nextflow Pipeline - Implemented Tool Commands](#ms_pipe)
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

<a name="rna_cs"></a>
### 1.3. Filtering RNA CS Control Sequences
**(1) Download the Sequences** 
Download the RNA CS sequences. Follow the link: [RNA CS (RNC)](https://nanoporetech.com/support/library-prep/RNA-and-cDNA/what-is-rna-cs-rcs)

**(2) Map the CS sequences against your raw sequencing data [minimap2](https://github.com/lh3/minimap2)**
```
minimap2 -ax map-ont rna_cs_control.fasta your_raw_reads.fasta > cs_mapped.sam
```

**(3) Extract read IDs**
First column in a SAM file holds the read ID and the 5th column holds the mapping quality (MAPQ). This command checks the first column that it is not a header row (starts with @) and checks each 5th row if it is greater than 0. From all rows complying to these conditions the read ID is written into a txt file.
For more info on SAM file format see here [`README_file_formats.md`](README_file_formats.md)
```
awk '$1 !~ /^@/ && $5 > 0 {print $1}' cs_mapped.sam > read_ids.txt
```

**(4) Exclude read IDs from the raw sequencing data with [seqkit](https://github.com/shenwei356/seqkit)**
Use seqkit to keep all raw reads whose read ID matches any in the read_id.txt file.
```
seqkit grep -v -f read_ids.txt raw_reads.fastq -o raw_reads_filtered.fastq
```

<a name="#pycoQC"></a>
### 1.4. pycoQC [https://a-slide.github.io/pycoQC/installation/]
The sequencing summary file to best use for this Quality Control report is availabel with the raq seqeuncing data. It is generated from the sequencing run. It hold information of ids, pore healt, start time, etc.
A report can also be generated from .fastq or .fasta files but it will show less information on the sequencing run.
```
pycoQC --summary file sequencing_summary.txt -o pycoQC_output.html 
```

<a name="phred_filter"></a>
### 1.5. chopper - Filtering based on PHRED-score [https://github.com/wdecoster/chopper]

```
gunzip –c xxx.fastq.gz | chopper --quality 10 --threads 18 --input | gzip xxx_filteredQ10.fastq.gz 

```

<a name="nanoComp"></a>
### 1.6. NanoComp [https://github.com/wdecoster/nanocomp]

```
# nanoComp report (compare fastq, fasta, summary)

NanoComp --summary xxx_summary.txt xxx_summary.txt xxx_summary.txt --names xxx_1 xxx_2 xxx_3 -o qc/nanocomp_xxx 
```

<a name="polya"></a>
### 1.7. Poly(A) tail trimming **[CUSTOM SCRIPT]**
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

seqkit grep -v -f read_ids.txt raw_reads.fastq -o raw_reads_noR.fastq
```

<a name="ms_pipe"></a>

<a name="nf_pipe"></a>
## 2. Nextflow Pipeline - Implemented Tool Commands 
## Transcriptome Reconstruction & Annotation

As the reconstruction and analyses of biological data can be very individual and to also provide the user with the knwoledge of how every step was used and in what order, all commands implemented into the nextflow pipeline are listed in this README. With this it is also easy to follow each step and/or use one ore more individually if necessary. 

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
...

<a name="busco"></a>
## 2.5  [BUSCO](https://busco.ezlab.org/busco_userguide.html) - Vertebrata Ortholog Assessment
```
# 4. BUSCO - Transcriptome Completeness Assessment
# offline:
busco -i xxx_transcripts.fa -l vertebrata_odb10 --download_path PATH/TO/busco_downloads/ -o OUTPUT_FILDER/ -m transcriptome --offline -c NO_THREADS
# with data download
busco -i xxx_transcripts.fa -l vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome -c NO_THREADS
```

<a name="transdec"></a>
## 2.6  [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) - ORF Prediction
```
# 5. TransDecoder - predict Open Reading Frames (ORFs)
TransDecoder.LongOrfs -t xxx_transcripts.fa 
```

<a name="eggnog"></a>
## 2.7  [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13) Annotation
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

