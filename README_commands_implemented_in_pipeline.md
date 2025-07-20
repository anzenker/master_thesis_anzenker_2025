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
    * 1.8. [ UniProt ](#up)  
+ 2. Fucntional Annotation Outside the Pipeline
    * 2.1. [ InterProScan ](#ips) 
+ 3. Custom Python Script for summarized GTF file of transcriptome and its annotations
    * 3.1. [ Final Output GTF](#fog)

   
###### ================================================================================

<a name="nf_pipe"></a>
## 1. Nextflow Pipeline - Implemented Tool Commands 
## Transcriptome Reconstruction & Annotation

As the reconstruction and analyses of biological data can be very individual and to also provide the user with the knwoledge of how every step was used and in what order, all commands implemented into the pipeline are listed in the following part. With this it is also easy to follow and use an individual step if necessary. 

#### **Commands implemented into the pipeline:**

<a name="mini_sam"></a>
### 1.1 minimap2 & samtools []
```
# 1. minimap2 - mapping of the raw reads to the reference genome
minimap2 -y --MD -ax splice -uf -k14 -t 24 xxx_genome.fa raw_reads.fastq.gz | samtools sort -@ 24 -o xxx_genome_mapped.bam

samtools index xxx_genome_mapped.bam
```

<a name="string2"></a>
### 1.2 stringtie2 []
```
# 2. stringtie2 - transcriptome reconstruction
stringtie -o stringtie2_xxx_transcripts.gtf -l NAME_PREFIX -L -p 24 xxx_genome_mapped.bam
```

<a name="gffread"></a>
### 1.3 gffread []
```
# 3. gffread - extract transcript sequences
gffread -w xxx_transcripts.fa -g xxx_genome.fa stringtie2_xxx_transcripts.gtf
```

<a name="busco"></a>
### 1.4 BUSCO []
```
# 4. BUSCO - Transcriptome Completeness Assessment
# offline:
busco -i xxx_transcripts.fa -l vertebrata_odb10 --download_path PATH/TO/busco_downloads/ -o OUTPUT_FILDER/ -m transcriptome --offline -c NO_THREADS
# with data download
busco -i xxx_transcripts.fa -l vertebrata_odb10 -o OUTPUT_FILDER/ -m transcriptome -c NO_THREADS
```

<a name="transdec"></a>
### 1.5 TransDecoder []
```
# 5. TransDecoder - predict Open Reading Frames (ORFs)
TransDecoder.LongOrfs -t xxx_transcripts.fa 
```

<a name="eggN"></a>
### 1.6 eggNOG []
EggNOG annotation ...
1. Download the required databases:
```
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
```
2. run the eggNOG command:
```
emapper.py  -m diamond --itype proteins -i xxx_transdecoder.pep -o NAME_PREFIX --data_dir /folder/with/databases
```

<a name="up"></a>
### 1.7 UniProt []
# 1.7 UniProt annotation **(planned)**


<a name="ips"></a>
## 2. Fucntional Annotation Outside the Pipeline
# 2.1. InterProScan annotation **(planned)**
It was decided to not include InterProScan into the Nextflow pipeline generated during this project. It caused memory issues while building the docker container for the pipeline. Additionally InterProScan provides an already available docker container which is easy to use. 
To install and run InterProScan with their docker container follow the commands listed below or visit the website [https://interproscan-docs.readthedocs.io/en/v5/HowToUseViaContainer.html] & the official docker website [https://hub.docker.com/r/interpro/interproscan].
```
docker pull interpro/interproscan:5.74-105.0
# ...
# Status: Downloaded newer image for interpro/interproscan:5.74-105.0
# docker.io/interpro/interproscan:5.74-105.0

# download InterProScan database (optional: md5sum check to verify the download)
curl -O http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/alt/interproscan-data-5.75-106.0.tar.gz.md5
# unpack download
tar -pxzf interproscan-data-5.75-106.0.tar.gz

# make working directories
mkdir input temp output

# copy fasta file
cp PATH/TO/transcriptome.fa /input

# command with explanations
#run InterProScan via docker
docker run --rm \
-u $(id -u):$(id -g) \                                          # run container as my current user, not root
-v $PWD/interproscan-5.75-106.0/data:/opt/interproscan/data \   # mount local folder into container
-v $PWD/work:/input \                                           # mount ... folder ...
-v $PWD/temp:/temp \                                            # mount ... folder ...
-v $PWD/output:/output interpro/interproscan:5.75-106.0 \       # mount ... folder ...
--input /input/xxx_canonical.fasta \                            # specify input file
--output-dir /output \
--tempdir /temp \
--cpu 18

# command to copy - adjust input file name
# if necessary adjust InterProScan version
docker run --rm -u $(id -u):$(id -g) -v $PWD/interproscan-5.75-106.0/data:/opt/interproscan/data -v $PWD/work:/input -v $PWD/temp:/temp -v $PWD/output:/output interpro/interproscan:5.75-106.0 --input /input/xxx.fasta --output-dir /output --tempdir /temp --cpu 18



```




<a name="fog"></a>
## 3. Custom Python Script for summarized GTF file of transcriptome and its annotations
## 3.1. Final output generation **(in progress)**

