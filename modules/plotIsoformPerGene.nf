#!/usr/bin/env nextflow

process plotIsoformPerGene {
    publishDir 'RESULTS/8_plots'
    container 'simp-test:latest'

    input:
    path gtf_input_file
    
    output:
    path "isoform_per_gene_barplot.png"

    script:
    """
    #!/opt/conda/bin/python

    import pandas as pd
    import matplotlib.pyplot as plt
    import re

    # col names in gtf file
    col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # read gtf into pd df
    df = pd.read_csv("$gtf_input_file", sep="\t", header=None, comment='#')
    df.columns = col_names
    # only keep rows with transcript feature
    transcript_df = df[df["feature"] == "transcript"].copy()


    # Extract gene_id and transcript_id using regular expressions
    def extract_id(text, key):
            match = re.search(f'{key} "([^"]+)"', text)
            return match.group(1) if match else None

    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))
    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))

    # Count how many transcripts (isoforms) each gene has
    isoform_counts = transcript_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    isoform_counts.columns = ["gene_id", "isoform_count"]

    # Count how many genes have 1, 2, 3… isoforms
    distribution = isoform_counts["isoform_count"].value_counts().sort_index()
    isoform_dist_df = distribution.reset_index()
    isoform_dist_df.columns = ["isoform_count", "number_of_genes"]

    # plot
    ax = isoform_dist_df.plot.bar(x='isoform_count', y='number_of_genes', rot=0, color="skyblue", legend=False)
    ax.set_xlabel("Number of Isoforms per Gene")
    ax.set_ylabel("Number of Genes")
    ax.set_title("Distribution of Isoform Count Per Genes")
    ax.bar_label(ax.containers[0], color="blue")
    plt.tight_layout()
    plt.savefig('isoform_per_gene_barplot.png')
    """
}