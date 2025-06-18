#  Commands Implemented into the pipeline

###### ================================================================================
###### OVERVIEW
+ 1. [Nextflow Pipeline - Implemented Tool Commands](#nf_pipe)
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

<a name="nf_pipe"></a>
## 2. Nextflow Pipeline Commands - Transcriptome Reconstruction & Annotation

As the reconstruction and analyses of biological data can be very individual and to also provide the user with the knwoledge of how every step was used and in what order, all commands implemented into the pipeline are listed in the following part. With this it is also easy to follow and use an individual step if necessary. 

#### Commands implemented into the pipeline:
```
# 1. minimap2 - mapping of the raw reads to the reference genome
minimap2 -y --MD -ax splice -uf -k14 -t 24 xxx_genome.fa raw_reads.fastq.gz | samtools sort -@ 24 -o xxx_genome_mapped.bam

samtools index xxx_genome_mapped.bam
```
```
# 2. stringtie2 - transcriptome reconstruction
stringtie -o stringtie2_xxx_transcripts.gtf -l NAME_PREFIX -L -p 24 xxx_genome_mapped.bam
```
```
# 3. gffread - extract transcript sequences
gffread -w xxx_transcripts.fa -g xxx_genome.fa stringtie2_xxx_transcripts.gtf
```
```
# 4. BUSCO - Transcriptome Completeness Assessment
# offline:
busco -i xxx_transcripts.fa -l /PATH/TO/busco_downloads/lineages/vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome --offline -c NO_THREADS
# with data download
busco -i xxx_transcripts.fa -l vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome -c NO_THREADS
```
```
# 5. TransDecoder - predict Open Reading Frames (ORFs)
TransDecoder.LongOrfs -t xxx_transcripts.fa 
```
```
# 6. EggNOG annotation
# download databases (~13G & ~9G & ~7G)
mkdir bin/
cd bin/

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
gunzip eggnog.db

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
gunzip eggnog_proteins.dmnd

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
tar –xvzf eggnog.taxa.tar.gz

cd ..

# run
emapper.py  -m diamond --itype proteins -i xxx_transdecoder.pep -o NAME_PREFIX --data_dir /folder/with/databases
```

# 7. InterProScan annotation **(planned)**
# 8. UniProt annotation **(planned)**
# 9. Final output generation **(in progress)**

