#!/usr/bin/env nextflow

process canonicalBestCov1 {

    publishDir 'RESULTS/6_canonical_transcriptome', mode:'copy'
    container 'simp-test:latest'

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
    publishDir 'RESULTS/6_canonical_transcriptome', mode:'copy'
    container 'simp-test:latest'

    input:
    path input_tsv_file
    path input_fasta_file

    output:
    path "${input_fasta_file.baseName}_canonical_ids.txt"

    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd
    df = pd.read_csv("$input_tsv_file", sep='\\t', header=None, names=["gene_id", "transcript_id", "coverage"])
    print(df.head())
    df["coverage"] = df["coverage"].astype(float)
    canonical_df = df.loc[df.groupby("gene_id")["coverage"].idxmax()]
    canonical_df["transcript_id"].to_csv("${input_fasta_file.baseName}_canonical_ids.txt", sep="\\t", header=None, index=False)
    """
}

process canonicalBestCov3 {

    publishDir 'RESULTS/6_canonical_transcriptome', mode:'copy'
    container 'simp-test:latest'

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