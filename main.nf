#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set default values for optional parameters - can be changed by user input 
params.help = false
params.outdir = 'RESULTS'
params.threads = 1
params.species_name = null
params.color = "#C79FEF" //lilac
params.skip_eggnog = false
params.skip_busco = false
params.skip_orf = false
params.no_plots = false


def helpMessage() {
        log.info """

        Usage:

        The typical command for running the pipeline is as follows:

        nextflow run ms_pipeline.nf --raw_reads xxx.fastq --genome xxx.fa --threads xx --with docker 


        Mandatory arguments:
            -r/--raw_reads              Path to input data - raw seqeuncing reads
            -g/--genome                 Path to input data - reference genome file (FASTA)
            -with-docker                run pipeline with docker

        Other options:
            -t/--threads                no. of threads (default: 1)
            -o/--outdir                 The output directory where the results will be saved
            -w/--work-dir               The temporary directory where intermediate data will be saved
            -sn/--species_name          For species names "A. arizonae", "A. marmoratus" and "A. neomexicanus" a specified plotting color is given.
            -c/--color                  Specific color for output plots created with Python. Default color is lavender (#C79FEF).

        
        Process Options:
            --skip_eggnog true          Skip EggNOG annotation step
            --skip_orf true             Skip TransDecoder ORF Prediction & eggNOG annotation.stripIndent
            --skip_busco true           Skip BUSCO analysis.
            --no_plots true             No output plot will be generated.

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
    stringtie --version     >> software_versions.txt 2>&1
    gffread --version       >> software_versions.txt 2>&1
    minimap2 --version      >> software_versions.txt 2>&1
    samtools --version      >> software_versions.txt 2>&1
    TransDecoder.LongOrfs   >> software_versions.txt 2>&1
    emapper.py --version    >> software_versions.txt 2>&1
    """
}

process buscoVertebrataCompleteness {
    publishDir "${params.outdir}/6_busco_vertebrata_completeness", mode: 'copy'

    input:
        val threads
        tuple path(input_fasta_file), val(label)
        path busco_downloads_path
    
    output:
        tuple path("busco_output_${input_fasta_file.baseName}/"), val(label)

    script:
    """
    busco -i $input_fasta_file -l vertebrata_odb10 --download_path $busco_downloads_path -o "busco_output_${input_fasta_file.baseName}" -m transcriptome --offline -c $threads  
    """
}

process canonicalBestCov1 {
    publishDir "${params.outdir}/4_canonical_transcriptome", mode:'copy'

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
    publishDir "${params.outdir}/4_canonical_transcriptome", mode:'copy'

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
    publishDir "${params.outdir}/4_canonical_transcriptome", mode:'copy'

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
    path "eggnog_${input_orf_pep.baseName}.emapper.hits", emit: hits
    path "eggnog_${input_orf_pep.baseName}.emapper.seed_orthologs", emit: seeds
    path "eggnog_${input_orf_pep.baseName}.emapper.annotations", emit: anno

    script:
    """
    emapper.py  -m diamond --itype proteins -i $input_orf_pep -o 'eggnog_${input_orf_pep.baseName}' --data_dir $eggDB_path --cpu $threads > emapper.log 2>&1

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
    when:
        !params.no_plots

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

process plotBUSCOCompleteness {
    when:
        !params.no_plots

    publishDir "${params.outdir}/8_plots"

    input:
    path python_script
    tuple path(full_table), val(label)
    val species_name

    output:
    path "busco_plot_${label}/6_busco_completeness_stacked_barplot.png"

    script:
    """
    #!/bin/bash
    
    python $python_script ${busco_full_table}/run_vertebrata_odb10/full_table.tsv $species_name $output_path busco_plot_${label}
    """
}

process plotIsoformPerGene {
    when:
        !params.no_plots

    publishDir "${params.outdir}/8_plots"

    input:
    path python_script
    path gtf_input_file
    val color
    
    output:
    path "2_ipg/${gtf_input_file.baseName}_isoform_per_gene_barplot.tsv"
    path "2_ipg/${gtf_input_file.baseName}_isoform_per_gene_barplot.png"
    path "2_ipg/${gtf_input_file.baseName}_isoform_per_gene_1_to_15_barplot.png"

    script:
    """
    #!/bin/bash
    
    python $python_script $gtf_input_file 2_ipg -plot_color "$color"
    """
}

process plotORFStatistics {
    when:
        !params.no_plots

    publishDir "${params.outdir}/8_plots", mode: 'copy'

    input:
    path python_script
    path input_fasta
    path input_pep
    val plot_color

    output:
    path "orf_compare_pivot.csv" 
    path "orf_compare_long.csv"
    path "orf_compare_grouped_bar.png"

    script:
    """
    #!/bin/bash

    python $python_script $input_fasta $input_pep -plot_color "$plot_color"
    """
}

process plotTotalTranscripts {
    when:
        !params.no_plots

    publishDir "${params.outdir}/8_plots", mode: 'copy'

    input:
    path python_script
    path input_fasta_1
    path input_fasta_2
    val plot_color
    
    output:
    path "transcript_counts_counts.csv"
    path "transcript_counts_grouped_bar.png"

    script:
    """
    #!/bin/bash
    
    python $python_script $input_fasta_1 $input_fasta_2 --color1 $plot_color
    """ 
}

process plotOverviewQuality {
    when:
        !params.no_plots

    publishDir "${params.outdir}/8_plots", mode: 'copy'

    input:
    path python_script
    path input_gtf
    path input_fasta
    path input_pep
    path input_busco
    path input_eggnog
    val species_name
    
    output:
    path "overview_quality/transcript_counts_counts.csv"
    path "overview_quality/transcript_counts_grouped_bar.png"

    script:
    """
    #!/bin/bash
    
    python $python_script $input_gtf $input_fasta $input_pep $input_busco $input_eggnog overview_quality --species_name $species_name
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
    publishDir "${params.outdir}/5_frame_selection", mode: 'copy'

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
    publishDir "${params.outdir}/NONE", mode: 'copy'

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

    //run processes
    //***************************************
    // 1. minimap & samtools
    // map raw reads to the genome assembly
    //***************************************
    minimap2RawToGenome(params.raw_reads, params.genome, params.threads)
    // --> % of mapping

    //***************************************
    // 2. stringtie
    // reconstruct transcriptome from mapping
    //***************************************
    stringtie2Transcriptome(params.threads, minimap2RawToGenome.out[0])
    if (!params.no_plots) {
        //----------------------------------------
        // plot isoform per gene
        //def isoform_plot_outdir = "${params.outdir}/8_plots"
        params.python_file_1 = '${projectDir}/python_scripts/2_plot_isoform_per_gene.py'
        //def python_script_path_2 = file('${projectDir}/python_scripts/2_plot_isoform_per_gene.py')
        plotIsoformPerGene(params.python_file_1, stringtie2Transcriptome.out, params.color)
        //----------------------------------------
    }
    //***************************************
    //3. gffread
    // read seqeuences with gtf and genome assembly into fasta file
    //***************************************
    gffreadToFasta(stringtie2Transcriptome.out, params.genome)
    // OUTPUT --> read metrics & plot read length distribution

    //***************************************
    // 4. awk, python, seqkit
    // generate canonical transcriptome from isoform transcripts with the highest coverage 
    //***************************************
    canonicalBestCov1(stringtie2Transcriptome.out) //.gtf --> .tsv
    canonicalBestCov2(canonicalBestCov1.out, gffreadToFasta.out) //.tsv, .fasta --> .txt
    canonicalBestCov3(canonicalBestCov2.out, gffreadToFasta.out) //.txt, .fasta --> .fasta

    // create channel from total transcriptome and canonical transcriptome
    transcriptome_ch = gffreadToFasta.out.combine(canonicalBestCov3.out).flatten()


    if (!params.no_plots) {
        //----------------------------------------
        // plot total transcriptome no all vs canonical 
        params.python_file_2 = '${projectDir}/python_scripts/4_plot_total_vs_canonical_transcript_count.py'
        //def python_script_path_1 = file('${projectDir}/python_scripts/4_plot_total_vs_canonical_transcript_count.py')
        plotTotalTranscripts(params.python_file_2, gffreadToFasta.out, canonicalBestCov3.out, params.color)
        //----------------------------------------
    }


    //transcriptome_ch = gffreadToFasta.out
    //                .combine(canonicalBestCov3.out)
    //                .map { a, b -> [a, b] } // ensure list structure
    //                .collectMany { it }    // flatten the list of pairs
    //                .enumerate()           // adds index


    // 5. TransDecoder (optional)
    if (!params.skip_orf) {
        //***************************************
        //5. ORF Prediction
        //***************************************
        transDecoderORF(canonicalBestCov3.out) //.fasta --> longest_orfs.pep
        //----------------------------------------
        if (!params.no_plots) {
            //----------------------------------------
            // plot ORF category distribution
            params.python_file_3 = '${projectDir}/python_scripts/5_plot_orf_statistics.py'
            //def python_script_path_3 = file('${projectDir}/python_scripts/5_plot_orf_statistics.py')
            plotORFStatistics(params.python_file_3, canonicalBestCov3.out, transDecoderORF.out, params.color)
            //----------------------------------------
        }   
    }

    // 6. BUSCO (otional)
    if (!params.skip_busco) {
        //***************************************
        //6. BUSCO Vertebrate - Completeness Assessment
        //***************************************
        params.busco_d_path = '${launchDir}/bin/busco_downloads/6_busco_completeness_stacked_barplot.py'
        //def busco_downloads_path = file('${launchDir}/bin/busco_downloads/6_busco_completeness_stacked_barplot.py')
        params.python_file_4 = '${projectDir}/python_scripts/6_busco_completeness_stacked_barplot.py'
        //def python_script_path_4 = file('${projectDir}/python_scripts/6_busco_completeness_stacked_barplot.py')

        buscoVertebrataCompleteness(params.threads, transcriptome_ch, params.busco_d_path)

        plotBUSCOCompleteness(params.python_file_4, buscoVertebrataCompleteness.out,   // tuple: (busco_full_table_<label>.tsv, label)
                                params.species_name ?: 'Species'
)
    }

    // 7. eggNOG (optional) - BUT only if ORF (5.) ran
        if (!params.skip_orf && !params.skip_eggnog) {
            //***************************************
            //7. EggNOG Annotation
            //***************************************
            params.eggnog_db = '${launchDir}/bin/'
            //def eggnog_databases_path = file('${launchDir}/bin/')
            eggnogAnnotation(transDecoderORF.out, params.threads, params.eggnog_db) //.pep, threads no., path --> .hits, .annotation, .seed_orthologs

        } 

    // 8) Overview (only if ALL present)
    if ( !params.no_plots && !params.skip_busco && !params.skip_orf && !params.skip_eggnog ) {
        params.python_file_5 = '${projectDir}/python_scripts/plot_pipeline_quality_overview_fucntioanlity.py'
        //def overview_py = file('${projectDir}/python_scripts/plot_pipeline_quality_overview_fucntioanlity.py')
        // choose which BUSCO channel to feed (e.g., canonical only):
        plotOverviewQuality(
        params.python_file_5,
        stringtie2Transcriptome.out,  // GTF
        gffreadToFasta.out,           // FASTA
        transDecoderORF.out,          // PEP
        buscoVertebrataCompleteness.out, // full_table.tsv
        eggnogAnnotation.out.hits,         // .annotations
        params.species_name ?: 'Species'
        )

    }
    //----------------------------------------
    workflow.onComplete {
        println "Pipeline finished successfully."
        println summary
    }

}
