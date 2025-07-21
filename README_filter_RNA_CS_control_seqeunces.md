#  Filter out RNA CS Control sequences from Raw Reads

###### ================================================================================
###### OVERVIEW
+ 1. [ Download the Sequences ](#down)
+ 2. [ Map the CS Sequences against your raw sequencing data ](#map)
+ 3. [ Extract read IDs ](#extract)
+ 4. [ Exclude read IDs from the raw sequencing reads ](#exclude)
###### ================================================================================

<a name="down"></a>
## 1. Download the Sequences 
Download the RNA CS sequences. Follow the link: [RNA CS (RNC)](https://nanoporetech.com/support/library-prep/RNA-and-cDNA/what-is-rna-cs-rcs)

<a name="map"></a>
## 2. Map the CS sequences against your raw sequencing data [minimap2](https://github.com/lh3/minimap2)
```
minimap2 -ax map-ont rna_cs_control.fasta your_raw_reads.fasta > cs_mapped.sam
```

<a name="extract"></a>
## 3. Extract read IDs
First column in a SAM file holds the read ID and the 5th column holds the mapping quality (MAPQ). This command checks the first column that it is not a header row (starts with @) and checks each 5th row if it is greater than 0. From all rows complying to these conditions the read ID is written into a txt file.
For more info on SAM file format see here [!README_file_formats.md]
```
awk '$1 !~ /^@/ && $5 > 0 {print $1}' cs_mapped.sam > read_ids.txt
```

<a name="exclude"></a>
## 4. Exclude read IDs from the raw sequencing data with [seqkit](https://github.com/shenwei356/seqkit)
Use seqkit to keep all raw reads whose read ID matches any in the read_id.txt file.
```
seqkit grep -v -f read_ids.txt raw_reads.fastq -o raw_reads_filtered.fastq
```

