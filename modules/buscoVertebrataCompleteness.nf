#!/usr/bin/env nextflow

process buscoVertebrataCompleteness {

    publishDir 'RESULTS/4_busco_vertebrata_completeness', mode: 'copy'
    container 'simp-test:latest'

    input:
        val threads
        path input_fasta_file
        path busco_downloads_path
    
    output:
        path "busco_output_${input_fasta_file.baseName}/"
        path "short_summary.${input_fasta_file.baseName}.vertebrata.odb12.XXX.txt"

    script:
    """
    busco -i $input_fasta_file -l vertebrata_odb10 --download_path $busco_downloads_path -o "busco_output_${input_fasta_file.baseName}" -m transcriptome --offline -c $threads

    cp busco_output_${input_fasta_file.baseName}/short_summary.*.txt short_summary.${input_fasta_file.baseName}.vertebrata.odb12.XXX.txt
    """
}