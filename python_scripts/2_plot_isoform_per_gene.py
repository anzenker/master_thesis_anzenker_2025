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

def extract_isoform_ids_per_gene(transcript_df, filename, output_path):
    
    # create a new column with the xtracted gene id from the attribute column
    transcript_df["gene_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "gene_id"))
    # create a new column with the xtracted transcript id from the attribute column
    transcript_df["transcript_id"] = transcript_df["attribute"].apply(lambda x: extract_id(x, "transcript_id"))

    # create a new df showing all transcript ids per gene id
    grouped_df = transcript_df.groupby(["seqname", "gene_id"])["transcript_id"].apply(list).reset_index()

    # add a count a isoform count column
    grouped_df["isoform_count"] = grouped_df["transcript_id"].apply(len)

    # Save to file
    output_file = f"{output_path}/{filename}_isoforms_per_gene_overview.tsv"
    grouped_df.to_csv(output_file, sep="\t", index=False)
    print(f"Isoform Counts saved to: {output_file}")


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
    isoform_dist_df.columns = ["isoform_count", "count_of_genes"]


    # save isoform per gene counts into txt file
    isoform_dist_df.to_csv(
        os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene.txt"),
        # the following option is used in the script used by the nextflow pipeline    
        #f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.txt",
        sep='\t', index= False
    )

    fig, ax = plt.subplots(figsize=(12, 6))

    # Ensure x is numeric and sorted
    isoform_dist_df["isoform_count"] = isoform_dist_df["isoform_count"].astype(int)
    isoform_dist_df = isoform_dist_df.sort_values("isoform_count")

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
    fig_1 = os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene_barplot.png")
    plt.savefig(fig_1)
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
    fig_2 = os.path.join(output_path, f"{os.path.splitext(filename)[0]}_isoform_per_gene_1_to_15_barplot.png")
    plt.savefig(fig_2)
    print(f"Plot saved to: {fig_2}")
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
    parser.add_argument("-plot_color", help="Provide color choice for plotting.", default="#C79FEF")
    args = parser.parse_args()

    # Extract basename from input file for naming of the final output plot 
    input_file_path = os.path.basename(args.input_gtf)
    filename = os.path.basename(os.path.realpath(args.input_gtf))

    # run functions
    transcripts_df = read_gtf_input_into_df(args.input_gtf)
    extract_isoform_ids_per_gene(transcripts_df, filename, args.output_path)
    plot_isoform_per_gene(transcripts_df, args.plot_color, filename, args.output_path)


    # the following option is used in the script used by the nextflow pipeline
    #plot_isoform_per_gene(transcripts_df, args.plot_color, filename)#, args.output_path)

if __name__ == "__main__":
    main()
