#!/usr/bin/env nextflow

process uniprotAnnotation {
    publishDir 'RESULTS/7_annotation', mode:'copy'
    container 'simp-test:latest'


    input:
    path input_orf_pep
    val threads
    path uniDB_path

    output:
    path "${input_orf_pep}_blastp.outfmt6"

    script:
    """
    blastp -query $input_orf_pep -db $uniDB_path  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $threads > '${input_orf_pep}_blastp.outfmt6'
    """
}