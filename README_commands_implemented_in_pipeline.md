#  Commands Implemented into the pipeline

###### ================================================================================
###### OVERVIEW
+ 1. [Nextflow Pipeline - Implemented Tool Commands](#nf_pipe)
    * 1.1. [ minimap2 & samtools ](#mini_sam)
    * 1.2. [ stringtie2 ](#string2)
    * 1.3. [ gffread ](#gffread)
    * 1.4. [ BUSCO ](#busco)
    * 1.5. [ TransDecoder ](#transdec)
    * 1.6. [ EggNOG ](#eggN)
    * 1.7. [ Custom Python Scripts ](#pyscr)
   
###### ================================================================================

<a name="nf_pipe"></a>
## 1. Nextflow Pipeline - Implemented Tool Commands 
## Transcriptome Reconstruction & Annotation

As the reconstruction and analyses of biological data can be very individual and to also provide the user with the knwoledge of how every step was used and in what order, all commands implemented into the nextflow pipeline are listed in this README. With this it is also easy to follow each step and/or use one ore more individually if necessary. 

#### **Commands implemented into the pipeline:**

<a name="mini_sam"></a>
### 1.1 [minimap2](https://lh3.github.io/minimap2/minimap2.html) & [samtools](https://www.htslib.org/doc/samtools.html)
```
# 1. minimap2 - mapping of the raw reads to the reference genome
minimap2 -y --MD -ax splice -uf -k14 -t 24 xxx_genome.fa raw_reads.fastq.gz | samtools sort -@ 24 -o xxx_genome_mapped.bam

samtools index xxx_genome_mapped.bam
```

<a name="string2"></a>
### 1.2 [stringtie2](https://github.com/skovaka/stringtie2)
```
# 2. stringtie2 - transcriptome reconstruction
stringtie -o stringtie2_xxx_transcripts.gtf -l NAME_PREFIX -L -p 24 xxx_genome_mapped.bam
```

<a name="gffread"></a>
### 1.3 [gffread](https://github.com/gpertea/gffread)
```
# 3. gffread - extract transcript sequences
gffread -w xxx_transcripts.fa -g xxx_genome.fa stringtie2_xxx_transcripts.gtf
```

<a name="busco"></a>
### 1.4 [BUSCO](https://busco.ezlab.org/busco_userguide.html)
```
# 4. BUSCO - Transcriptome Completeness Assessment
# offline:
busco -i xxx_transcripts.fa -l vertebrata_odb10 --download_path PATH/TO/busco_downloads/ -o OUTPUT_FILDER/ -m transcriptome --offline -c NO_THREADS
# with data download
busco -i xxx_transcripts.fa -l vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome -c NO_THREADS
```

<a name="transdec"></a>
### 1.5 [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)
```
# 5. TransDecoder - predict Open Reading Frames (ORFs)
TransDecoder.LongOrfs -t xxx_transcripts.fa 
```

<a name="eggN"></a>
### 1.6 [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13)
1. Download the database
For database download instruction refer to [README_commands_used.md](README_commands_used.md)

2. run the eggNOG command:
```
emapper.py  -m diamond --itype proteins -i xxx_transdecoder.pep -o NAME_PREFIX --data_dir /folder/with/databases
```

<a name="pyscr"></a>
### 1.7 Custom Python Scripts
- plot Isoform per Gene
- plot ORF distribution
- plot BUSCO Completeness
- plot eggNOG annotation %
