# Examplary Results - Pipeline Output per Step

###### ================================================================================
###### OVERVIEW
+ 1. [ minimap2RawtoGenome ](#one)
+ 2. [ stringtie2Transcriptome ](#two)
+ 3. [ gffreadToFasta ](#three)
+ 4. [ canonicalBestCov ](#four)
+ 5. [ transDecoderORF ](#five)
+ 6. [ buscoVertebrataCompleteness ](#six)
+ 7. [ eggnogAnnotation ](#seven)
+ 8. [ Plots ](#eight)
    + 8.1 [ plotIsoformPerGene ](#nine)
    + 8.2 [ plotTotalTranscripts ](#ten)
    + 8.3 [ plotORFStatistics ](#eleven)
    + 8.4 [ plotBUSCOCompleteness ](#twelve) 
    + 8.5 [ plotOverviewQuality ](#thirteen)

================================================================================

<a name="one"></a>
### 1. minimap2RawtoGenome
- BAM file
- indexed BAM.BAI file

<a name="two"></a>
### 2. stringtie2Transcriptome
- GTF file of the total Transcriptome

<a name="three"></a>
### 3. gffreadToFasta
- FASTA file of the total Transcriptome

<a name="four"></a>
### 4. canonicalBestCov
- transcript_ids_and_covergae.tsv
- canonical_transcript_ids.txt
- FASTA file with canonical transcripts

<a name="five"></a>
### 5. transDecoderORF
Only done with the canonical Transcriptome.
- PEP file with the longest predicted ORFs

<a name="six"></a>
### 6. buscoVertebrataCompleteness
- busco output for total Transcriptome
busco output for canonical Transcriptome

<a name="seven"></a>
### 7. eggnogAnnotation
- eggNOG annotation for canonical Transcriptome

<a name="eigth"></a>
### 8. Plots
<a name="nine"></a>
### 8.1 plotIsoformPerGene

![2_isoform_per_gene_barplot.png](/images/2_isoform_per_gene_barplot.png)

![2_isoform_per_gene_1_to_15_barplot.png](/images/2_isoform_per_gene_1_to_15_barplot.png)

<a name="ten"></a>
### 8.2 plotTotalTranscripts

![4_plot_total_vs_canonical_counts.png](/images/4_plot_total_vs_canonical_counts.png)

<a name="eleven"></a>
### 8.3 plotORFStatistics

![5_plot_orf_statistics.png](/images/5_plot_orf_statistics.png)

<a name="twelve"></a>
### 8.4 plotBUSCOCompleteness

![total Transcriptome: 6_busco_completeness_stacked_barplot.png](/images/6_busco_completeness_stacked_barplot_total.png)

![canonical Transcriptome: 6_busco_completeness_stacked_barplot.png](/images/6_busco_completeness_stacked_barplot_canonical.png)

<a name="thirteen"></a>
### 8.5 plotOverviewQuality

![pipeline_transcriptome_quality_overview_functionality.png](/images/pipeline_transcriptome_quality_overview_functionality.png)
