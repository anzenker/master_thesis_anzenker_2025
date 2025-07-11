#!/usr/bin/env nextflow

process stringtie2Transcriptome {

    publishDir 'RESULTS/2_stringtie2_transcriptome'
    container 'simp-test:latest'

    input:
        val threads
        path input_bam_file

    output:
        path "stringtie2_${input_bam_file.baseName}_transcripts.gtf"
    
    script:
    """
    stringtie -o "stringtie2_${input_bam_file.baseName}_transcripts.gtf" -l "${input_bam_file.baseName}" -L -p $threads $input_bam_file
    """
}