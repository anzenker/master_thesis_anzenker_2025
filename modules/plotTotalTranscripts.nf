#!/usr/bin/env nextflow

process plotTotalTranscripts {
    publishDir 'RESULTS/8_plots'
    container 'simp-test:latest'

    input:
    path python_script
    path input_fasta_file
    path input_pep_file
    
    output:
    path "no_all_transcripts_vs_canonical_transcripts.png"

    script:
    """
    #!/bin/bash
    
    python $python_script $input_fasta_file $input_pep_file
    """ 
}