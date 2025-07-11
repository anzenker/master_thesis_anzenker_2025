#!/usr/bin/env nextflow

process minimap2RawToGenome {

    publishDir 'RESULTS/1_minimap2_output', mode: 'copy'
    container 'simp-test:latest'

    input:
        path raw_reads
        path genome_file
        val threads

    output:
        path "${raw_reads.baseName}_mapped.bam"
        path "${raw_reads.baseName}_mapped.bam.bai"

    script:
    """
    minimap2 -y --MD -ax splice -uf -k14 -t $threads $genome_file $raw_reads | samtools sort -@ $threads -o "${raw_reads.baseName}_mapped.bam"

    samtools index "${raw_reads.baseName}_mapped.bam" -o "${raw_reads.baseName}_mapped.bam.bai"
    """
}