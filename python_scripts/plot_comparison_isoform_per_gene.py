import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse
from datetime import datetime
import os

def read_gtf_input_into_df(gtf_input):
    # col names in gtf file
    col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # read gtf into pd df
    df = pd.read_csv(gtf_input, sep="\t", header=None, comment='#')
    df.columns = col_names
    # only keep rows with transcript feature
    transcript_df = df[df["feature"] == "transcript"].copy()

    # Extract gene_id and transcript_id
    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))
    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))
    
    return transcript_df

# Extract gene_id and transcript_id using regular expressions
def extract_id(text, key):
        match = re.search(f'{key} "([^"]+)"', text)
        return match.group(1) if match else None

def calculate_isoform_distribution(transcript_df, label):
    isoform_counts = transcript_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    isoform_counts.columns = ["gene_id", "isoform_count"]
    distribution = isoform_counts["isoform_count"].value_counts().sort_index()
    dist_df = distribution.reset_index()
    dist_df.columns = ["isoform_count", label]
    return dist_df

def plot_isoform_per_gene(merged, output_plot):

    # Plot
    colors = ["#8e6c88", "#b99666", "#1469A7"]  #neo, marm, ari
    ax = merged.plot(
        x="isoform_count",
        kind="bar",
        stacked=False,
        color=colors,
        width=0.75
    )
    ax.set_xlabel("Number of Isoforms per Gene")
    ax.set_ylabel("Number of Genes")
    ax.set_title("Isoform Count Distribution per Species")
    ax.legend(title="Species")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot isoform per gene from GTF output file from stringtie")
    parser.add_argument("input_gtf_neo", help="Input GTF file - A. neomexicanus")
    parser.add_argument("input_gtf_marm", help="Input GTF file - A. marmoratus")
    parser.add_argument("input_gtf_ari", help="Input GTF file - A. arizonae")   
    parser.add_argument("output_plot", help="Output plot showing isoforms per gene - Comparison between 3 species")

    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{os.path.splitext(args.output_plot)[0]}_log_{timestamp}.txt"

    # dataframes
    transcripts_df_neo = read_gtf_input_into_df(args.input_gtf_neo)
    transcripts_df_marm = read_gtf_input_into_df(args.input_gtf_marm)
    transcripts_df_ari = read_gtf_input_into_df(args.input_gtf_ari)

    # Calculate distributions
    dist_neo = calculate_isoform_distribution(transcripts_df_neo, "neomexicanus")
    dist_marm = calculate_isoform_distribution(transcripts_df_marm, "marmoratus")
    dist_ari = calculate_isoform_distribution(transcripts_df_ari, "arizonae")

    # Merge
    merged = pd.merge(dist_neo, dist_marm, on="isoform_count", how="outer")
    merged = pd.merge(merged, dist_ari, on="isoform_count", how="outer")
    merged = merged.fillna(0).sort_values(by="isoform_count")
    merged["isoform_count"] = merged["isoform_count"].astype(int)


    plot_isoform_per_gene(merged, args.output_plot)

if __name__ == "__main__":
    main()
