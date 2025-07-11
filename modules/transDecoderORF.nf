#!/usr/bin/env nextflow

//also include TransDecoder.Predict with homology options????
process transDecoderORF {

    publishDir 'RESULTS/5_frame_selection', mode:'copy'
    container 'simp-test:latest'

    input:
    path input_fasta_file

    output:
    //path "base_freqs.dat"
    //path "__checkpoints_longorfs"
    //path "longest_orfs.cds"
    //path "longest_orfs.gff3"
    path "${input_fasta_file}.transdecoder_dir/longest_orfs.pep"

    script:
    """
    TransDecoder.LongOrfs -t $input_fasta_file
    """
}