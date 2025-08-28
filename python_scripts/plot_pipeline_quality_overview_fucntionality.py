
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import textwrap
import os

def read_eggnog_into_df(eggnog_input_path): 
    # read file into DataFrame
    eggNOG_df = pd.read_csv(eggnog_input_path, sep="\t", header=None)

    # only keep first two columns (transcript ID & Gene ID)
    eggnog_ids_df = pd.concat([eggNOG_df[0], eggNOG_df[1]], axis=1)
    eggnog_ids_df.columns = ['qseqid', 'sseqid_eggNOG']
    annotated_transcript = set(eggnog_ids_df['qseqid'])
    # extract just the transcript ID part before '.px' - make compatible also with BUSCO
    #eggnog_ids_df['qseqid'] = eggnog_ids_df["qseqid"].str.replace(r'\.p\d+$', '', regex=True)

    return annotated_transcript

def read_busco_into_df(busco_input_path): 
    # read file into DataFrame
    col_names = ['sseqid_BUSCO', 'Status', 'qseqid', 'Score', 'Length', 'OrthoDB', 'url', 'Description']
    busco_df = pd.read_csv(busco_input_path, sep='\t', comment='#', names=col_names)
    #busco_df = busco_df.dropna()
    #busco_df.columns = col_names

    busco_counts = busco_df['Status'].value_counts()

    count_single_complete = len(busco_df[busco_df['Status'] == 'Complete'])
    count_duplicated_complete = len(busco_df[busco_df['Status'] == 'Duplicated'])
    count_complete_all = count_single_complete + count_duplicated_complete
    count_total = len(busco_df)

    p_complete = round((count_complete_all / count_total) * 100, 2)

    return busco_df, p_complete

def parse_transdecoder_pep(input_pep):

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

    transcripts_with_orf = set(max_orfs.keys())
    total_with_orf = len(transcripts_with_orf)

    return total_with_orf



def pipeline_plot_transcriptome_quality_overview(species_name, species_color, count_canonical, count_orf_predicted, count_eggnog_anno, percent_busco_complete, output_path):
    # Labels and counts
    labels = ['Canonical', 'With ORF', 'With \n eggNOG \n Annotation']
    counts = [count_canonical, count_orf_predicted, count_eggnog_anno]

    # Percentages relative to canonical
    percentages = [
        100,
        round((count_orf_predicted / count_canonical) * 100, 2),
        round((count_eggnog_anno / count_canonical) * 100, 2)
    ]

    fig, ax1 = plt.subplots(figsize=(8, 6))

    # Plot counts on primary y-axis
    bars = ax1.bar(labels, counts, color=species_color, width=0.4)
    ax1.set_ylabel("Count Transcripts", fontsize=12)
    ax1.set_ylim(0, max(counts) * 1.2)

    # Add count labels
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2, height + max(counts) * 0.01,
                f"{str(height)} \n {round((height/count_canonical)*100, 2)}%", ha='center', va='bottom', fontsize=10)

    # Secondary y-axis for percentage
    ax2 = ax1.twinx()
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_ylabel("Percentage of BUSCOs found in Transcriptome", fontsize=12, color='gray')
    ax2.set_yticks([y / 100 * max(counts) for y in range(0, 110, 10)])
    ax2.set_yticklabels([f"{y}%" for y in range(0, 110, 10)], color='gray')

    # Add BUSCO line
    busco_height = percent_busco_complete / 100 * max(counts)
    ax1.bar('Complete \n BUSCOs', busco_height , width=0.4, color='gray')
    ax1.text(len(labels), busco_height + max(counts) * 0.01,
            f"Percent \n BUSCO \n Complete: \n {percent_busco_complete}%",
            ha='center', va='bottom', fontsize=10, color='gray')

    plt.title(fr"Transcriptome Quality Overview - $\it{{{species_name}}}$", fontsize=14)
    plt.tight_layout()

    # Save
    fig = os.path.join(output_path, "pipeline_transcriptome_quality_overview_functionality.png")
    plt.savefig(fig)
    print(f"Plot saved to: {fig}")
    
    plt.close()




def default_color_for(species: str) -> str:
    # override here for your 3 common species
    if species.strip() == "A. marmoratus":  return "#b99666"
    if species.strip() == "A. arizonae":    return "#1469A7"
    if species.strip() == "A. neomexicanus":return "#688e26"
    return "#C79FEF"  # fallback

def main():

    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input file
    parser.add_argument("input_gtf", help="Path to stringtie GTF file")
    parser.add_argument("input_canonical_ids", help="Path to TXT file with IDs from Canonical Transcripts")
    parser.add_argument("input_pep", help="Path to PEP file with Open Reading Frame Prediction")
    parser.add_argument("input_busco", help="Path to BUSCO full_table.tsv")
    parser.add_argument("input_eggnog", help="Path to eggNOG.annotation file")

    parser.add_argument("output_path", help="Path to save output plot image to.")
    parser.add_argument("--species_name", help="Species Name")
    args = parser.parse_args()

    if args.species_name:
        species_name = args.species_name
    else:
        species = "Species"

    col_names = ["BUSCO_id", "Status", "Sequence", "Score", "Length", "OrthoDB", "url", "Description"]


    ids_all = pd.read_csv(args.input_gtf)
    count_transc_all = len(ids_all)

    ids_canonical = pd.read_csv(args.input_canonical_ids)
    count_canonical = len(ids_canonical)

    count_with_orf = parse_transdecoder_pep(args.input_pep)


    eggnog_ids = read_eggnog_into_df(args.input_eggnog)
    busco_ids, p_complete = read_busco_into_df(args.input_busco)

    pipeline_plot_transcriptome_quality_overview(species_name,
                                        default_color_for(species_name),
                                        count_canonical, 
                                        count_with_orf, 
                                        len(eggnog_ids), 
                                        p_complete, args.output_path
        )
if __name__ == "__main__":
    main()

