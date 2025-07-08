from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import textwrap

def plot_orf_distribution(input_fasta, input_pep):

    total_no_transcripts = len([1 for line in open(input_fasta) if line.startswith(">")])

    orf_dict = defaultdict(list)
    orf_type_counter = Counter()

    # Parse .pep file and group by transcript base ID
    for record in SeqIO.parse(input_pep, "fasta"):
        desc_parts = record.description.split()
        orf_type = "unknown"
        for part in desc_parts:
            if part.startswith("type:"):
                orf_type = part.split(":")[1]
                break
        transcript_id = record.id.rsplit(".", 1)[0]  # AC73_xxx.p1/p2 → AC73_xxx
        seq_length = len(record.seq)
        # save for each transcript id the each detected type of ORF with their length
        orf_dict[transcript_id].append((orf_type, seq_length))

    # only choose one ORF per transcript id --> complete & longest
    max_orfs = {}

    for transcript_id, orfs in orf_dict.items():
        # Choose the longest ORF regardless of type
        best_orf = max(orfs, key=lambda x: x[1])  # x[1] is the length
        max_orfs[transcript_id] = best_orf
        orf_type_counter[best_orf[0]] += 1

    # Compute transcript counts
    transcripts_with_orf = set(max_orfs.keys())
    total_with_orf = len(transcripts_with_orf)
    single_orf = sum(1 for orfs in orf_dict.values() if len(orfs) == 1)
    multi_orf = sum(1 for orfs in orf_dict.values() if len(orfs) > 1)

    # Adjust based on your transcriptome
    total_transcripts = total_no_transcripts  # ← update this to actual total transcript number
    total_with_orf = len(max_orfs)
    transcripts_without_orf = total_transcripts - total_with_orf

    # Prepare counts and labels
    categories = [
        "All Transcripts",
        "With ORF",
        "Without ORF",
        "1 ORF",
        ">1 ORF",
        "Complete ORF",
        "5' Partial ORF",
        "3' Partial ORF",
        "Internal ORF"
    ]
    counts = [
        total_transcripts,
        total_with_orf,
        transcripts_without_orf,
        single_orf,
        multi_orf,
        orf_type_counter.get("complete", 0),
        orf_type_counter.get("5prime_partial", 0),
        orf_type_counter.get("3prime_partial", 0),
        orf_type_counter.get("internal", 0)
    ]
    # Define a color for each category (you can customize these)
    colors = [
        "blue",  # All Transcripts
        "blue",  # With ORF
        "blue",  # Without ORF
        "lightblue",  # 1 ORF
        "lightblue",  # >1 ORF
        "green",  # Complete ORF
        "orange",  # 5' Partial ORF
        "tan",  # 3' Partial ORF
        "grey"   # Internal ORF
    ]

    percentages = [f"{(c/total_transcripts)*100:.1f}%" for c in counts]
    numbers = [f"{total_transcripts}" for c in counts]
    # Combine number and percentage labels
    labels = [f"{count}\n({percent})" for count, percent in zip(counts, percentages)]

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(categories, counts, color=colors)

    # Add labels above each bar
    for bar, label in zip(bars, labels):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.,
            height,
            label,
            ha='center',
            va='bottom',
            fontsize=10
        )

    # Axis labels and formatting
    ax.set_ylabel("Number of Transcripts / ORFs")
    ax.set_title("Transcript and ORF Statistics")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("orf_prediction.png")

if __name__ == "__main__":
    # script description
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input files given by user 
    parser.add_argument("input_fasta" , help="Path to Transcriptome FASTA file")
    parser.add_argument("input_pep", help="Path to TransDecoder PEP file")
    args = parser.parse_args()

    # run plotting
    plot_orf_distribution(args.input_fasta, args.input_pep)
