#!/usr/bin/env nextflow

process gffreadToFasta {

    publishDir 'RESULTS/3_gffread_transcriptome'
    container 'simp-test:latest'

    input:
        path input_gtf_file
        path genome_file

    output:
        path "${input_gtf_file.baseName}.fasta"
    
    script:
    """
    gffread -w "${input_gtf_file.baseName}.fasta" -g $genome_file $input_gtf_file
    """

}