#!/usr/bin/env nextflow

// set default values for optional parameters - can be changed by user input 
params.help = false
params.outdir = 'RESULTS'
params.threads = 1
params.skip_eggnog = false
params.help = false

//include { minimap2RawToGenome } from './modules/minimap2RawToGenome.nf'
//include { stringtie2Transcriptome } from './modules/stringtie2Transcriptome.nf'
//include { gffreadToFasta } from './modules/gffreadToFasta.nf'
//include { canonicalBestCov1 } from './modules/canonicalBestCov.nf'
//include { canonicalBestCov2 } from './modules/canonicalBestCov.nf'
//include { canonicalBestCov3 } from './modules/canonicalBestCov.nf'
//include { buscoVertebrataCompleteness } from './modules/buscoVertebrataCompleteness.nf'
//include { transDecoderORF } from './modules/transDecoderORF.nf'
//include { eggnogAnnotation } from './modules/eggnogAnnotation.nf'
//include { uniprotAnnotation } from './modules/uniprotAnnotation.nf'

//include { plotIsoformPerGene } from './modules/plotIsoformPerGene.nf'
//include { plotLengthDistribution } from './modules/plotLengthDistribution.nf'
//include { plotTotalTranscripts } from './modules/plotTotalTranscripts.nf'

def helpMessage() {
        log.info """

        Usage:

        The typical command for running the pipeline is as follows:

        nextflow run ms_pipeline.nf --raw_reads xxx.fastq --genome xxx.fa --threads xx --with docker 


        Mandatory arguments:
            --raw_reads                 Path to input data - raw seqeuncing reads
            --genome                    Path to input data - reference genome file (FASTA)
            --with-docker               run pipeline with docker

        Other options:
            --threads                   no. of threads
            --outdir                    The output directory where the results will be saved
            -w/--work-dir               The temporary directory where intermediate data will be saved
        
        Process Options:

            --skip_eggnog               Skip EggNOG annotation step
            --skip_blast                Skip blastp annotation step
        """.stripIndent()
    }



/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions.txt'

    script:
    """
    conda --version         > software_versions.txt
    busco --version         >> software_versions.txt 2>&1
    seqkit --version        >> software_versions.txt 2>&1
    python --version        >> software_versions.txt 2>&1
    pip --version           >> software_versions.txt 2>&1
    biopython --version     >> software_versions.txt 2>&1
    pandas --version        >> software_versions.txt 2>&1
    bbmap --version         >> software_versions.txt 2>&1
    blast --version         >> software_versions.txt 2>&1
    augustus --version      >> software_versions.txt 2>&1
    metaeuk --version       >> software_versions.txt 2>&1
    prodigal --version      >> software_versions.txt 2>&1
    hmmer --version         >> software_versions.txt 2>&1
    sepp --version          >> software_versions.txt 2>&1
    r-base --version        >> software_versions.txt 2>&1         
    r-ggplot2 --version     >> software_versions.txt 2>&1
    eggnog-mapper --version >> software_versions.txt 2>&1
    """
}

process buscoVertebrataCompleteness {
    publishDir "${params.outdir}/7_busco_vertebrata_completeness", mode: 'copy'

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

process canonicalBestCov1 {
    publishDir "${params.outdir}/5_canonical_transcriptome", mode:'copy'

    input: 
    path input_gtf_file

    output:
    path "transcript_ids_and_coverage.tsv"

    script:
    """
    awk ' BEGIN { OFS="\t" } { if (\$3 == "transcript") {print \$10, \$12, \$14} }' $input_gtf_file | sed 's/[\";"]//g' > transcript_ids_and_coverage.tsv
    """
}

process canonicalBestCov2 {
    publishDir "${params.outdir}/5_canonical_transcriptome", mode:'copy'

    input:
    path input_tsv_file
    path input_fasta_file

    output:
    path "${input_fasta_file.baseName}_canonical_ids.txt"

    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd

    # read input file
    df = pd.read_csv("$input_tsv_file", sep='\\t', header=None, names=["gene_id", "transcript_id", "coverage"])
    
    # verify value type
    df["coverage"] = df["coverage"].astype(float)

    # keep isoform of transcripts per gene which has the max covergae value
    canonical_df = df.loc[df.groupby("gene_id")["coverage"].idxmax()]

    # save transcript ids into txt file
    canonical_df["transcript_id"].to_csv("${input_fasta_file.baseName}_canonical_ids.txt", sep="\\t", header=None, index=False)
    """
}

process canonicalBestCov3 {
    publishDir "${params.outdir}/5_canonical_transcriptome", mode:'copy'

    input:
    path input_tsv_file
    path input_fasta_file

    output:
    path "${input_fasta_file.baseName}_canonical.fasta"

    script:
    """
    seqkit grep -f $input_tsv_file $input_fasta_file -o "${input_fasta_file.baseName}_canonical.fasta"
    """
}

process eggnogAnnotation {
    //when:
    //!params.skip_eggnog

    publishDir "${params.outdir}/7_eggnog_annotation", mode:'copy'

    input:
    path input_orf_pep
    val threads
    path eggDB_path

    output:
    path "eggnog_${input_orf_pep.baseName}.emapper.hits"
    path "eggnog_${input_orf_pep.baseName}.emapper.seed_orthologs"
    path "eggnog_${input_orf_pep.baseName}.emapper.annotations"

    script:
    """
    emapper.py  -m diamond --itype proteins -i $input_orf_pep -o 'eggnog_${input_orf_pep.baseName}' --data_dir $eggDB_path --cpu $threads

    """
}

process gffreadToFasta {

    publishDir "${params.outdir}/3_gffread_transcriptome"

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

process minimap2RawToGenome {

    publishDir "${params.outdir}/1_minimap2_output", mode: 'copy'

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

process plotIsoformPerGene {
    publishDir "${params.outdir}/8_plots"

    input:
    path python_script
    path gtf_input_file
    
    output:
    path "${gtf_input_file.baseName}_isoform_per_gene_barplot.png"
    path "${gtf_input_file.baseName}_isoform_per_gene_barplot.txt"


    script:
    """
    #!/bin/bash
    
    python $python_script $gtf_input_file
    """
}
process plotLengthDistribution {
    publishDir "${params.outdir}/8_plots", mode: 'copy'

    input:
    path input_fasta_file
    
    output:
    path "${input_fasta_file.baseName}_length_distribution.png"

    script:
    """
    #!/opt/conda/bin/python
    
    from Bio import SeqIO
    import pandas as pd 
    import matplotlib.pyplot as plt

    records = list(SeqIO.parse("$input_fasta_file", "fasta"))
    lengths = [len(record.seq) for record in records]
    no_seqs = len(records)
    
    df = pd.DataFrame({"Sequence Length": lengths})
    length_value_counts_df = df.value_counts().reset_index()
    length_value_counts_df.columns = ["Transcript Length", "Count"]

    ax = length_value_counts_df.plot.bar(x='Count', y='Transcript Length', rot=0, color="skyblue")
    ax.set_xlabel("Number of Transcripts")
    ax.set_ylabel("Transcript Length")
    ax.set_title("Length Distribution of Transcriptome ($input_fasta_file)")
    plt.tight_layout()
    plt.savefig('${input_fasta_file.baseName}_length_distribution.png')
    """ 
}

process plotTotalTranscripts {
    publishDir "${params.outdir}/8_plots", mode: 'copy'

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

process stringtie2Transcriptome {
    
    publishDir "${params.outdir}/2_stringtie2_transcriptome", mode: 'copy'

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

//also include TransDecoder.Predict with homology options????
process transDecoderORF {
    publishDir "${params.outdir}/6_frame_selection", mode: 'copy'

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

// does not work due to memory issues
process uniprotAnnotation {
    publishDir "${params.outdir}/5_frame_selection", mode: 'copy'
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

workflow {

    // show help message and exit
    if (params.help){
        helpMessage()
        exit 0
    }
    // parameter check
    if (!params.raw_reads) {
        exit 1, "Missing paramter: --raw_reads"
    }
    if (!params.genome){
        exit 1, "Missing parameter: --genome"
    }

    // check input file paths 
    Channel.fromPath(params.raw_reads, checkIfExists: true)
            .ifEmpty {
                exit 1, "File ${params.raw_reads} cannot be found.}"
            }
    Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty {
                exit 1, "File ${params.genome} cannot be found."
            }

    // pipeline run info
    log.info "==========PIPELINE START=========="
    def summary = [:]
    summary['Pipeline Name']    = 'ms_pipeline'
    summary['Run Name']         = workflow.runName
    summary['Raw Reads']        = params.raw_reads
    summary['Reference Genome'] = params.genome
    summary['Threads']          = params.threads
    summary['Output Directory'] = params.outdir
    summary['User']             = workflow.userName
    summary['Start Time']       = workflow.start

    // define parameters & input paths
    def eggnog_databases_path = file('./bin/')
    def uniprot_databases_path = file('./bin/uniprot_sprot.fasta')
    def busco_downloads_path = file('./bin/busco_downloads/')
    def python_script_path_1 = file('./python_scripts/plot_transcript_numbers.py')
    def python_script_path_2 = file('./python_scripts/plot_isoform_per_gene.py')
    def python_script_path_3 = file('./python_scripts/plot_transcript_numbers.py')


    //run processes

    //***************************************
    //1. map raw reads to the genome assembly
    //***************************************
    minimap2RawToGenome(params.raw_reads, params.genome, params.threads)
    // --> % of mapping

    //***************************************
    //2. reconstruct transcriptome from mapping
    //***************************************
    stringtie2Transcriptome(params.threads, minimap2RawToGenome.out[0])
    // --> plot isoform per gene

    //***************************************
    //3. read seqeuences with gtf and genome assembly into fasta file
    //***************************************
    gffreadToFasta(stringtie2Transcriptome.out, params.genome)
    // OUTPUT --> read metrics & plot read length distribution



    //***************************************
    //4. generate canonical transcriptome from isoform transcripts with the highest coverage 
    //***************************************
    canonicalBestCov1(stringtie2Transcriptome.out) //.gtf --> .tsv
    canonicalBestCov2(canonicalBestCov1.out, gffreadToFasta.out) //.tsv, .fasta --> .txt
    canonicalBestCov3(canonicalBestCov2.out, gffreadToFasta.out) //.txt, .fasta --> .fasta
    // OUTPUT --> read metrics & plot read length distribution
    // --> ?????overlaying read length distribution??????

    // create channel from total transcriptome and canonical transcriptome
    transcriptome_ch = gffreadToFasta.out
                        .combine(canonicalBestCov3.out)
                        .flatten()
                        .eachWithIndex( { fasta, idx -> tuple(fasta, idx) } )

    //***************************************
    //5. ORF Prediction
    //***************************************
    transDecoderORF(canonicalBestCov3.out) //.fasta --> longest_orfs.pep


    //***************************************
    //6. BUSCO Vertebrate - Completeness Assessment
    //***************************************
    buscoVertebrataCompleteness(params.threads, transcriptome_ch, busco_downloads_path)


    //***************************************
    //7. EggNOG Annotation
    //***************************************
    eggnogAnnotation(transDecoderORF.out, params.threads, eggnog_databases_path) //.pep, threads no., path --> .hits, .annotation, .seed_orthologs

    //***************************************
    //8. PLOTS
    //***************************************
    //----------------------------------------
    // plot isoform per gene
    plotIsoformPerGene(python_script_path_2, stringtie2Transcriptome.out)
    //----------------------------------------
    //----------------------------------------
    // plot transcriptome length distribution - combine .fasta into channel and run process
    plotLengthDistribution(transcriptome_ch)
    //----------------------------------------
    //----------------------------------------
    // plot comparison of total transcript count
    plotTotalTranscripts(python_script_path_1, transDecoderORF.out, canonicalBestCov3.out)
    //----------------------------------------
    //----------------------------------------
    // plot BUSCO output

    //----------------------------------------    
    //----------------------------------------
    // plot ORF category distribution

    //----------------------------------------
    //----------------------------------------
    // plot eggNOG annotation results overview
    
    //----------------------------------------

}