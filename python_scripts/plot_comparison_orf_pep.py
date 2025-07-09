from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import textwrap
import pandas as pd

def count_total_number_in_fasta(input_fasta):
    return sum(1 for line in open(input_fasta) if line.startswith(">"))

def group_pep_by_transcript(input_pep, total_transcripts, label):

    orf_dict = defaultdict(list)
    orf_type_counter = Counter()

    for record in SeqIO.parse(input_pep, "fasta"):
        desc_parts = record.description.split()
        orf_type = "unknown"
        for part in desc_parts:
            if part.startswith("type:"):
                orf_type = part.split(":")[1]
                break
        transcript_id = record.id.rsplit(".", 1)[0]
        seq_length = len(record.seq)
        orf_dict[transcript_id].append((orf_type, seq_length))

    max_orfs = {}
    for transcript_id, orfs in orf_dict.items():
        best_orf = max(orfs, key=lambda x: x[1])
        max_orfs[transcript_id] = best_orf
        orf_type_counter[best_orf[0]] += 1

    total_with_orf = len(max_orfs)
    total_without_orf = total_transcripts - total_with_orf
    single_orf = sum(1 for orfs in orf_dict.values() if len(orfs) == 1)
    multi_orf = sum(1 for orfs in orf_dict.values() if len(orfs) > 1)

    complete_orf = orf_type_counter.get("complete", 0)
    fivePartial_orf = orf_type_counter.get("5prime_partial", 0)
    threePartial_orf = orf_type_counter.get("3prime_partial", 0)
    internal_orf = orf_type_counter.get("internal", 0)

    counts_row = {
        "Species": label,
        "Total Transcripts": total_transcripts,
        "With ORF": total_with_orf,
        "Without ORF": total_without_orf,
        "1 ORF": single_orf,
        ">1 ORF": multi_orf,
        "Complete ORF": complete_orf,
        "5' Partial ORF": fivePartial_orf,
        "3' Partial ORF": threePartial_orf,
        "Internal ORF": internal_orf
    }

    return pd.DataFrame([counts_row])

def bar_plot_orf(df):

    df["Species"] = pd.Categorical(df["Species"], 
                                categories=["neomexicanus", "marmoratus", "arizonae"],
                                ordered=True)

    df_long = df.melt(id_vars="Species", var_name="Category", value_name="Count")

    # Define desired order of x-axis categories
    category_order = [
        "Total Transcripts",
        "With ORF",
        "Without ORF",
        "1 ORF",
        ">1 ORF",
        "Complete ORF",
        "5' Partial ORF",
        "3' Partial ORF",
        "Internal ORF"
    ]

    # Apply order to the Category column
    df_long["Category"] = pd.Categorical(df_long["Category"], categories=category_order, ordered=True)

    pivot_df = df_long.pivot(index="Category", columns="Species", values="Count")

    species_colors = {
        "neomexicanus": "#688e26",
        "marmoratus": "#b99666",
        "arizonae": "#1469A7"
    }
    colors = [species_colors[species] for species in pivot_df.columns]

    ax = pivot_df.plot(kind="bar", color=colors, width=0.95, figsize=(12, 6))

    plt.ylabel("Number of Transcripts / ORFs")
    plt.title("ORF Categories per Species")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("orf_categories_grouped_bar.png")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        Plot ORF category distributions for multiple species.'''),
        usage='python %(prog)s <args>', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("input_fasta_neo", help="Transcriptome FASTA - A. neomexicanus")
    parser.add_argument("input_pep_neo", help="TransDecoder PEP - A. neomexicanus")
    parser.add_argument("input_fasta_marm", help="Transcriptome FASTA - A. marmoratus")
    parser.add_argument("input_pep_marm", help="TransDecoder PEP - A. marmoratus")
    parser.add_argument("input_fasta_ari", help="Transcriptome FASTA - A. arizonae")
    parser.add_argument("input_pep_ari", help="TransDecoder PEP - A. arizonae")
    args = parser.parse_args()

    total_neo = count_total_number_in_fasta(args.input_fasta_neo)
    total_marm = count_total_number_in_fasta(args.input_fasta_marm)
    total_ari = count_total_number_in_fasta(args.input_fasta_ari)

    df_neo = group_pep_by_transcript(args.input_pep_neo, total_neo, "neomexicanus")
    df_marm = group_pep_by_transcript(args.input_pep_marm, total_marm, "marmoratus")
    df_ari = group_pep_by_transcript(args.input_pep_ari, total_ari, "arizonae")

    full_df = pd.concat([df_neo, df_marm, df_ari], ignore_index=True)
    bar_plot_orf(full_df)