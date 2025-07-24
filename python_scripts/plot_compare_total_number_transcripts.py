import matplotlib.pyplot as plt
import argparse
import textwrap

def plot_transcript_numbers(input_fasta, input_fasta_canonical, plot_color):
    #count total number of transcripts in FASTA files
    total_no_transcripts = len([1 for line in open(input_fasta) if line.startswith(">")])
    total_no_canonical = len([1 for line in open(input_fasta_canonical) if line.startswith(">")])

    # plotting
    categories = ["All Transcripts", "Canonical Transcripts"]
    counts = [total_no_transcripts, total_no_canonical]
    percentages = [f"{(c / total_no_transcripts) * 100:.1f}%" for c in counts]
    labels = [f"{count}\n({percent})" for count, percent in zip(counts, percentages)]
    colors = [plot_color] * 2

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    bars = ax.bar(categories, counts, color=colors)

    # Add labels above each bar
    for bar, label in zip(bars, labels):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height, label, ha='center', va='bottom', fontsize=10
        )

    # set further plot info
    ax.set_ylabel("Number of Transcripts")
    ax.set_title("All Transcripts vs Canonical Transcripts")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig("no_all_transcripts_vs_canonical_transcripts.png")

if __name__ == "__main__":
    # script description
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input files given by user 
    parser.add_argument("input_fasta" , help="Path to Transcriptome FASTA file")
    parser.add_argument("input_fasta_canonical", help="Path to Canonical Transcriptome FASTA file")
    parser.add_argument("-plot_color", required=False, help="Provide color choice for plotting.", default="#C79FEF")

    args = parser.parse_args()

    plot_color = args.plot_color if args.plot_color is not None else "#C79FEF"
    # run plotting
    plot_transcript_numbers(args.input_fasta, args.input_fasta_canonical, plot_color)
