#!/usr/bin/env nextflow

process eggnogAnnotation {
    publishDir 'RESULTS/7_annotation', mode:'copy'
    container 'simp-test:latest'


    input:
    path input_orf_pep
    val threads
    path eggDB_path

    output:
    path "eggnog_${input_orf_pep.baseName}.emapper.hits"
    path "eggnog_${input_orf_pep.baseName}.emapper.seed_orthologs"
    path "eggnog_${input_orf_pep.baseName}.annotations"

    script:
    """
    emapper.py  -m diamond --itype proteins -i $input_orf_pep -o 'eggnog_${input_orf_pep.baseName}' --data_dir $eggDB_path --cpu $threads

    """
}