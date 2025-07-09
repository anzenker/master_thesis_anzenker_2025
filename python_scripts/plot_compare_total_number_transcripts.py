import matplotlib.pyplot as plt
import argparse
import textwrap
import gzip
import pandas as pd
import seaborn as sns

def open_func(input_fasta):
    gzip.open if input_fasta.endswith(".gz") else open

def count_total_number_in_fasta(input_fasta, input_canonical_fasta, label):
    total_transcripts = len([1 for line in open(input_fasta) if line.startswith(">")])
    total_canonical = len([1 for line in open(input_canonical_fasta) if line.startswith(">")])

    data = [ 
        ["Total Transcripts", total_transcripts, label], 
        ["Total Canonical", total_canonical, label]
    ]
    df = pd.DataFrame(data)
    df.columns = ["Total", "Count", "Species"]
    
    return df

def plot_transcript_numbers(df):

    
    df["Species"] = pd.Categorical(
        df["Species"],
        categories=["neomexicanus", "marmoratus", "arizonae"],
        ordered=True
    )

    df["Total"] = pd.Categorical(
        df["Total"],
        categories=["Total Transcripts", "Total Canonical"],
        ordered=True
    )

    # Pivot the data so we get transcript types as index and species as columns
    pivot_df = df.pivot(index="Total", columns="Species", values="Count")

    # Define custom colors per species (must match the column order!)
    species_colors = {
        "neomexicanus": "#688e26",
        "marmoratus": "#b99666",
        "arizonae": "#1469A7"
    }
    colors = [species_colors[species] for species in pivot_df.columns]

    # Plot
    ax = pivot_df.plot(
        kind="bar",
        width=0.75,
        color=colors,
        figsize=(8, 6)
    )

    # Customization
    ax.set_ylabel("Number of Transcripts")
    ax.set_title("Total vs Canonical Transcripts per Species")
    ax.set_xlabel("")
    plt.xticks(rotation=0)
    plt.legend(title="Species")
    plt.tight_layout()
    plt.savefig("grouped_transcripts_by_species.png")
    plt.show()

if __name__ == "__main__":
    # script description
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input files given by user 
    parser.add_argument("input_fasta_neo" , help="Path to Transcriptome FASTA file - A. neomexicanus")
    parser.add_argument("input_fasta_canonical_neo", help="Path to Canonical Transcriptome FASTA file - A. neomexicanus")
    parser.add_argument("input_fasta_marm" , help="Path to Transcriptome FASTA file - A. marmoratus")
    parser.add_argument("input_fasta_canonical_marm", help="Path to Canonical Transcriptome FASTA file - A. marmoratus")
    parser.add_argument("input_fasta_ari" , help="Path to Transcriptome FASTA file - A. arizonae")
    parser.add_argument("input_fasta_canonical_ari", help="Path to Canonical Transcriptome FASTA file - A. arizonae")
    args = parser.parse_args()

    # calc the total numbers
    df_neo = count_total_number_in_fasta(args.input_fasta_neo, args.input_fasta_canonical_neo, "neomexicanus")
    df_marm = count_total_number_in_fasta(args.input_fasta_marm, args.input_fasta_canonical_marm, "marmoratus")
    df_ari = count_total_number_in_fasta(args.input_fasta_ari, args.input_fasta_canonical_ari, "arizonae")

    #merge
    # concatenate all rows into one long dataframe
    merged = pd.concat([df_neo, df_marm, df_ari], ignore_index=True)

    # run plotting
    plot_transcript_numbers(merged)
