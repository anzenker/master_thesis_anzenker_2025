#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.outdir = 'RESULTS'


include { sayHello } from './modules/sayHello.nf'
include { minimap2RawToGenome } from './modules/minimap2RawToGenome.nf'
include { stringtie2Transcriptome } from './modules/stringtie2Transcriptome.nf'
include { gffreadToFasta } from './modules/gffreadToFasta.nf'
include { canonicalBestCov1 } from './modules/canonicalBestCov.nf'
include { canonicalBestCov2 } from './modules/canonicalBestCov.nf'
include { canonicalBestCov3 } from './modules/canonicalBestCov.nf'
include { buscoVertebrataCompleteness } from './modules/buscoVertebrataCompleteness.nf'
include { transDecoderORF } from './modules/transDecoderORF.nf'
include { eggnogAnnotation } from './modules/eggnogAnnotation.nf'
include { uniprotAnnotation } from './modules/uniprotAnnotation.nf'

include { plotIsoformPerGene } from './modules/plotIsoformPerGene.nf'
include { plotLengthDistribution } from './modules/plotLengthDistribution.nf'
include { plotTotalTranscripts } from './modules/plotTotalTranscripts.nf'

// DOES NOT WORK LIKE THIS
process plot_busco_completeness {

    publishDir '${outdir_path}/plots', mode: 'copy'
    container 'simp-test:latest'

    input:
        path summary_txt_files

    output:
        path 'busco_completeness_all_transcriptomes_plot.png'

    script:
    """
    mkdir summaries
    cp $summary_txt_files summaries/
    python /Users/annikazenker/Desktop/NEXTFLOW_tests/python_scripts/generate_plot.py -wd summaries
    """
}


workflow {
    def eggnog_databases_path = file('./bin/')
    def uniprot_databases_path = file('./bin/uniprot_sprot.fasta')
    def busco_downloads_path = file('./bin/busco_downloads/')
    def python_script_path_1 = file('./python_scripts/plot_transcript_numbers.py')
    params.greeting = "Bonjour!"


    //create a channel for inputs
    greeting_ch = Channel.of(params.greeting, 'Bonjour', 'Konnichiwa')
    greetings_array = ["Tschau", "Holá"]
    greeting_ch_2 = Channel.of(greetings_array).flatten()
    //run processes
    //sayHello(greeting_ch)
    sayHello(greeting_ch_2)
    
    //1. map raw reads to the genome assembly
    minimap2RawToGenome(params.raw_reads, params.genome_file, params.threads)

    //2. reconstruct transcriptome from mapping
    stringtie2Transcriptome(params.threads, minimap2RawToGenome.out[0])
    //3. read seqeuences with gtf and genome assembly into fasta file
    gffreadToFasta(stringtie2Transcriptome.out, params.genome_file)
    //4. generate canonical transcriptome from isoform transcripts with the highest coverage 
    canonicalBestCov1(stringtie2Transcriptome.out) //.gtf --> .tsv
    canonicalBestCov2(canonicalBestCov1.out, gffreadToFasta.out) //.tsv, .fasta --> .txt
    canonicalBestCov3(canonicalBestCov2.out, gffreadToFasta.out) //.txt, .fasta --> .fasta


    transcriptome_ch = gffreadToFasta.out
                            .combine(canonicalBestCov3.out)
                            .flatten()
                            .eachWithIndex( { fasta, idx -> tuple(fasta, idx) } ).view()

    transcriptome_ch.view()




    plotIsoformPerGene(stringtie2Transcriptome.out)
    plotLengthDistribution(transcriptome_ch)

    buscoVertebrataCompleteness(params.threads, transcriptome_ch, busco_downloads_path)

    //busco_vertebrata_completeness.out[1].view().collect().view().flatten().view()

    // Pass to plot process
    //plot_busco_completeness(busco_vertebrata_completeness.out[1].collect().flatten())
    // plotting script not working bc of problem with R dependency!!!!

    // collects the inputs into an array
    transDecoderORF(canonicalBestCov3.out.collect())

    //eggnogAnnotation(transDecoderORF.out, params.threads, eggnog_databases_path )

    //uniprotAnnotation(transDecoderORF.out, params.threads, uniprot_databases_path )

    plotTotalTranscripts(python_script_path_1, transDecoderORF.out, canonicalBestCov3.out)

}