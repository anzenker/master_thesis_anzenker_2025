import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse
from datetime import datetime
import os

def read_gtf_input_into_df(gtf_input):
    # col names in gtf file
    col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # read gtf into pandas DataFrame
    df = pd.read_csv(gtf_input, sep="\t", header=None, comment='#')
    df.columns = col_names
    # only keep rows with transcript feature
    transcript_df = df[df["feature"] == "transcript"].copy()

    return transcript_df

def extract_id(text, key):

    match = re.search(f'{key} "([^"]+)"', text)
    return match.group(1) if match else None

def extract_isoform_ids_per_gene(transcript_df, filename, output_path):
    
    # create a new column with the extracted gene id from the attribute column
    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))
    # create a new column with the extracted transcript id from the attribute column
    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))

    # create a new df showing all transcript ids per gene id
    grouped_df = transcript_df.groupby(["seqname", "gene_id"])["transcript_id"].apply(list).reset_index()

    # add a isoform count column
    grouped_df["isoform_count"] = grouped_df["transcript_id"].apply(len)

    # save
    tab_1 = os.path.join(output_path, "2_isoforms_per_gene.tsv")
    grouped_df.to_csv(tab_1, sep="\t", index=False)
    print(f"Isoform Counts saved to: {tab_1}")


def plot_isoform_per_gene(transcript_df, plot_color, filename, output_path):
# the following option is used in the script used by the nextflow pipeline
#def plot_isoform_per_gene(transcript_df, plot_color, filename):#, output_path):

    # extract the gene ids and transcript ids from the 'attributes' column 
    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))

    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))

    # Count how many transcripts (isoforms) each gene has
    isoform_counts = transcript_df.groupby("gene_id")["transcript_id"].nunique().reset_index()

    isoform_counts.columns = ["gene_id", "isoform_count"]

    # count how many genes have 1, 2, 3… isoforms
    distribution = isoform_counts["isoform_count"].value_counts().sort_index()

    isoform_dist_df = distribution.reset_index()
    isoform_dist_df.columns = ["isoform_count", "count_of_genes"]


    # save
    tab_2 = os.path.join(output_path, "4_isoforms_per_gene.txt")
    isoform_dist_df.to_csv(tab_2, sep="\t", index=False)
    print(f"Isoform Counts saved to: {tab_2}")

    fig, ax = plt.subplots(figsize=(12, 6))

    isoform_dist_df["isoform_count"] = isoform_dist_df["isoform_count"].astype(int)
    isoform_dist_df = isoform_dist_df.sort_values("isoform_count")

    #save
    isoform_dist_df.to_csv(os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.tsv"))

    x = isoform_dist_df["isoform_count"]
    y = isoform_dist_df["count_of_genes"]

    bars = ax.bar(x, y, color=plot_color, width=0.8)

    # Fill missing isoform_count values with 0 to ensure all values appear
    min_count = isoform_dist_df["isoform_count"].min()
    max_count = isoform_dist_df["isoform_count"].max()
    #ax.set_xticks(x)  # set x-ticks to actual numeric isoform counts
    ax.set_xticks(range(min_count, max_count + 1, 2))  # show every 2nd tick
    ax.set_xlabel("Count of Isoforms per Gene")
    ax.set_ylabel("Count of Genes")
    ax.set_title("Distribution of Isoform Count Per Genes", fontsize=12)
    plt.suptitle(f"input:file {filename}", fontsize=10)
    ax.bar_label(bars, color='black', fontsize=8)
    plt.tight_layout()

    # Save
    fig_1 = os.path.join(output_path, f"2_isoform_per_gene_barplot.png")
    plt.savefig(fig_1, dpi=500)
    print(f"Plot saved to: {fig_1}")
    #plt.close()
    
    bars = ax.bar(x, y, color=plot_color, width=0.8)

    ax.set_xticks(range(0, 16, 1))  # show every 2nd tick
    plt.xlim(0, 15.5)  # or ax.set_xlim(1, 15)
    ax.set_xlabel("Count of Isoforms per Gene")
    ax.set_ylabel("Count of Genes")
    ax.set_title("Distribution of Isoform Count Per Genes \n zoomed in on range 1 to 15", fontsize=12)
    plt.suptitle(f"input:file {filename}", fontsize=10)
    ax.bar_label(bars, color='black', fontsize=8)
    plt.tight_layout()

    # Save
    fig_2 = os.path.join(output_path, f"2_isoform_per_gene_1_to_15_barplot.png")
    plt.savefig(fig_2)
    print(f"Plot saved to: {fig_2}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot isoform per gene from GTF output file from stringtie")
    parser.add_argument("input_gtf", help="Input GTF file")

    parser.add_argument("output_path", help="Path to save output plot image to.")
    parser.add_argument("-plot_color", help="Provide color choice for plotting.", default="#C79FEF")
    args = parser.parse_args()

    # Extract basename from input file for naming of the final output plot 
    input_file_path = os.path.basename(args.input_gtf)
    filename = os.path.basename(os.path.realpath(args.input_gtf))

    # run functions
    transcripts_df = read_gtf_input_into_df(args.input_gtf)
    extract_isoform_ids_per_gene(transcripts_df, filename, args.output_path)
    plot_isoform_per_gene(transcripts_df, args.plot_color, filename, args.output_path)

if __name__ == "__main__":
    main()
