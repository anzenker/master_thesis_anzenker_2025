import pandas as pd
import matplotlib.pyplot as plt
import re
import argparse
from datetime import datetime
import os

def read_gtf_input_into_df(gtf_input):
    """
    Reads a GTF file and returns pandas DataFrame containing all entries of the feature 'transcript'.
    
    Parameters
    ----------
    gtf_input : str or path 

    Returns
    ----------
    transcript_df : pandas.DataFrame
    """
    # col names in gtf file
    col_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    # read gtf into pandas DataFrame
    df = pd.read_csv(gtf_input, sep="\t", header=None, comment='#')
    df.columns = col_names
    # only keep rows with transcript feature
    transcript_df = df[df["feature"] == "transcript"].copy()

    return transcript_df

def extract_id(text, key):
    """
    Extracts the value of a specified key from a GTF attribute string using regular expression.

    Paramters
    ----------
    text: str 
        A string representing the GTF attributes field.
    key: str 
        Key to extract the value for, 'gene_id' or 'transcript_id'

    Returns:
    str or None 
        Value associated with the key or None if the key is not found.
    ----------

    """
    match = re.search(f'{key} "([^"]+)"', text)
    return match.group(1) if match else None

def plot_isoform_per_gene(transcript_df, plot_color, filename, output_path):
# the following option is used in the script used by the nextflow pipeline
#def plot_isoform_per_gene(transcript_df, plot_color, filename):#, output_path):
    """
    Plots the isoform count per gene ID and saves it as a bar plot.

    Parameters
    ----------
    transcript_df : pandas.DataFrame
        DataFrame holding GTF transcript entries with an 'attributes' column.
    plot_color : str
        Specified color for the plot.
    basename : str
        Basename of input file for naming of the output plot image.
    output_path : str
        Path to save the output plot image to.

    Return
    ----------
    None
        Saves the isoform counts per gene into a tab separated txt file format.
        Saves the output plot image.

    """
    # Extract the gene ids and transcript ids from the 'attributes' column 
    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))
    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))

    # Count how many transcripts (isoforms) each gene has
    isoform_counts = transcript_df.groupby("gene_id")["transcript_id"].nunique().reset_index()
    isoform_counts.columns = ["gene_id", "isoform_count"]

    # Count how many genes have 1, 2, 3… isoforms
    distribution = isoform_counts["isoform_count"].value_counts().sort_index()
    isoform_dist_df = distribution.reset_index()
    isoform_dist_df.columns = ["isoform_count", "number_of_genes"]

    # save isoform per gene counts into txt file
    isoform_dist_df.to_csv(
        os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene.txt"),
        # the following option is used in the script used by the nextflow pipeline    
        #f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.txt",
        sep='\t', index= False
    )

    # plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax = isoform_dist_df.plot.bar(x='isoform_count', y='number_of_genes', rot=0, color=plot_color, legend=False)
    ax.set_xlabel("Number of Isoforms per Gene")
    ax.set_ylabel("Number of Genes")
    ax.set_title("Distribution of Isoform Count Per Genes", fontsize=12)
    plt.suptitle(f"input:file {filename}", fontsize=10)
    ax.bar_label(ax.containers[0], color='black', fontsize=8)
    plt.tight_layout()
    #Save Plot
    plt.savefig(os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.png"))
    # the following option is used in the script used by the nextflow pipeline
    #plt.savefig(f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.png")
    plt.close()

def main():
    """
    Main function to take user specific input and run the isoform per gene plotting script. 

    This function:
    - Parses arguments for input GTF file, output path, and plot color.
    - Extracts the filename.
    - Loads GTF data into a DataFrame.
    - Generates and saves a a txt file and a bar plot of isoform counts per gene.
    """
    parser = argparse.ArgumentParser(description="Plot isoform per gene from GTF output file from stringtie")
    parser.add_argument("input_gtf", help="Input GTF file")
    # commented out in the nextflow pipeline
    parser.add_argument("output_path", help="Path to save output plot image to.")
    parser.add_argument("plot_color", help="Provide color choice for plotting.")
    args = parser.parse_args()

    # Extract basename from input file for naming of the final output plot 
    input_file_path = os.path.basename(args.input_gtf)
    filename = os.path.basename(os.path.realpath(args.input_gtf))

    # run functions
    transcripts_df = read_gtf_input_into_df(args.input_gtf)
    plot_isoform_per_gene(transcripts_df, args.plot_color, filename, args.output_path)
    # the following option is used in the script used by the nextflow pipeline
    #plot_isoform_per_gene(transcripts_df, args.plot_color, filename)#, args.output_path)

if __name__ == "__main__":
    main()
