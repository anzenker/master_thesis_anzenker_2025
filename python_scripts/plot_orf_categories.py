from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import textwrap
from matplotlib import colors as mcolors

def parse_input_files(input_fasta, input_pep):
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

    return total_no_transcripts, best_orf, max_orfs, orf_type_counter, orf_dict

def plot_orf_distribution(total_no_transcripts, best_orf, max_orfs, orf_dict, plot_color):
    # Compute transcript counts
    transcripts_with_orf = set(max_orfs.keys())
    total_with_orf = len(transcripts_with_orf)
    single_orf = sum(1 for orfs in orf_dict.values() if len(orfs) == 1)
    multi_orf = sum(1 for orfs in orf_dict.values() if len(orfs) > 1)

    # Adjust based on your transcriptome
    total_transcripts = total_no_transcripts  # ← update this to actual total transcript number
    total_with_orf = len(max_orfs)
    transcripts_without_orf = total_transcripts - total_with_orf

    categories = [
        "Count Transcriptome (T)",
        "Count T. with ORF",
        "Count T. without ORF",
        "Count T. with 1 ORF",
        "Count T. with > 1 ORF"
    ]

    counts = [
        total_transcripts,
        total_with_orf,
        transcripts_without_orf,
        single_orf,
        multi_orf,
    ]

    colors = [plot_color] * 5

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
    plt.savefig("orf_prediction_categories.png")

def plot_pie_chart_orf_categories(orf_type_counter, plot_color):

    # Prepare counts and labels
    categories = [
        "Complete ORF",
        "5' Partial ORF",
        "3' Partial ORF",
        "Internal ORF"
    ]
    counts = [
        orf_type_counter.get("complete", 0),
        orf_type_counter.get("5prime_partial", 0),
        orf_type_counter.get("3prime_partial", 0),
        orf_type_counter.get("internal", 0)
    ]

    colors = generate_shades(plot_color, 4)


    fig, ax = plt.subplots()
    ax.pie(counts, labels=categories, colors=colors, autopct='%1.1f%%', startangle=180, labeldistance=1.1, pctdistance=.9)
    ax.set_title("ORF Type Distribution")
    plt.tight_layout()
    plt.show()
    plt.savefig("pie_chart_orf_categories.png")
    plt.close()

def generate_shades(base_color, n):
    base_rgb = mcolors.to_rgb(base_color)
    shades = [tuple(min(1, c + i*0.15) for c in base_rgb) for i in range(n)]
    return shades

if __name__ == "__main__":
    # script description
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input files given by user 
    parser.add_argument("-input_fasta", help="Path to Transcriptome FASTA file")
    parser.add_argument("-input_pep", help="Path to TransDecoder PEP file")
    parser.add_argument("-plot_color", default="#C79FEF", help="Color for plotting in hex code, default = lilac (#C79FEF)")
    args = parser.parse_args()

    plot_color = args.plot_color if args.plot_color is not None else  "#C79FEF"
    # run plotting
    total_no_transcripts, best_orf, max_orfs, orf_type_counter, orf_dict = parse_input_files(args.input_fasta, args.input_pep)

    #best_orf not yet used 

    plot_orf_distribution(total_no_transcripts, best_orf, max_orfs, orf_dict, plot_color)
    plot_pie_chart_orf_categories(orf_type_counter, plot_color)
