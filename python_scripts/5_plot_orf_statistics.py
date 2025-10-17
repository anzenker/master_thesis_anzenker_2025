from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import argparse, textwrap, os
import pandas as pd

def parse_one_species(fasta_path: str, pep_path: str):
    # returns orf categories (complete, 5', 3', internal) from TransDecoder
    total_transcripts = sum(1 for line in open(fasta_path) if line.startswith(">"))

    # collect all ORFs per transcript base id
    orf_dict = defaultdict(list)
    for rec in SeqIO.parse(pep_path, "fasta"):
        # TransDecoder description has 'type:<complete|5prime_partial|3prime_partial|internal>'
        orf_type = "unknown"
        for part in rec.description.split():
            if part.startswith("type:"):
                orf_type = part.split(":", 1)[1]
                break
        # remove .p1 suffix (AC73_xxx.p1 -> AC73_xxx)
        base_tid = rec.id.rsplit(".", 1)[0]  
        orf_dict[base_tid].append((orf_type, len(rec.seq)))

    # keep the longest ORF per transcript & count types
    longest_orf_per_tx = {}
    orf_type_counter = Counter()
    for tid, orfs in orf_dict.items():
        # longest
        best = max(orfs, key=lambda x: x[1])    
        longest_orf_per_tx[tid] = best
        orf_type_counter[best[0]] += 1

    with_orf_count = len(longest_orf_per_tx)
    return total_transcripts, with_orf_count, orf_type_counter

def summarize_species_row(species: str, fasta: str, pep: str):
    total, with_orf, c = parse_one_species(fasta, pep)
    row = {
        "species": species,
        "with_orf": with_orf,
        "complete":        c.get("complete", 0),
        "5partial":        c.get("5prime_partial", 0),
        "3partial":        c.get("3prime_partial", 0),
        "internal":        c.get("internal", 0),
        "total_transcripts": total
    }
    return row

def plot_grouped_from_pivot(pivot_df: pd.DataFrame, colors_map: dict, output_path: str):
    # pivot table --> index = orf category, columns = species, values = counts

    # order categories
    order = ["total_transcripts", "with_orf", "complete", "5partial", "3partial", "internal"]
    order = [c for c in order if c in pivot_df.index]
    pivot_df = pivot_df.loc[order]

    # colors in column order
    bar_colors = [colors_map.get(sp, "#888888") for sp in pivot_df.columns]

    ax = pivot_df.T.plot(kind="bar", color=None, figsize=(8, 5))  # bars grouped by category for each species
    # we want categories on x-axis; above transposed, so revert:
    plt.clf()
    ax = pivot_df.plot(kind="bar", figsize=(9,5), width=0.5, color=bar_colors, edgecolor='white')

    ax.set_ylabel("Count")
    ax.set_xlabel("Categories")
    ax.set_title("ORF categories per species")
    ax.legend(title="Species", frameon=False)


    totals = pivot_df.loc["total_transcripts"]  # Series: per-species totals

    ymax = pivot_df.max().max() * 1.10
    ax.set_ylim(0, ymax)

    for col_idx, container in enumerate(ax.containers):
        denom = float(totals.iloc[col_idx]) if len(totals) > col_idx and totals.iloc[col_idx] else 0.0
        for bar in container:
            h = float(bar.get_height())
            if h > 0:
                pct = (h / denom * 100.0) if denom > 0 else 0.0
                ax.text(
                    bar.get_x() + bar.get_width()/2,
                    h + ymax*0.01,
                    f"{pct:.1f}%\n{int(h)}",
                    ha="center", va="bottom", fontsize=9
                )
                
    plt.xticks(rotation=0)
    plt.tight_layout()

    # Save
    fig = os.path.join(output_path, "5_plot_orf_statistics.png")
    plt.savefig(fig, dpi=500)
    print(f"Plot saved to: {fig}")

    plt.close()

def default_color_for(species: str) -> str:
    if species.strip() == "A. marmoratus":  return "#b99666"
    if species.strip() == "A. arizonae":    return "#1469A7"
    if species.strip() == "A. neomexicanus":return "#688e26"
    return "#b99666"  

def main():

    jobs = []
    species_order = []
    color_map = {}

    #set colors per species
    if args.plot_color1 is not None:
        col1 = args.plot_color1
    elif args.species_name1 == "A. marmoratus":
        col1 = "#b99666"
    elif args.species_name1 == "A. arizonae":
        col1 = "#1469A7"      
    elif args.species_name1 == "A. neomexicanus":
        col1 = "#688e26"
    else:
        col1 = "#b99666"
    
    if args.input_fasta2 is not None:
        if args.plot_color2 is not None:
            col2 = args.plot_color2
        elif args.species_name2 == "A. marmoratus":
            col2 = "#b99666"
        elif args.species_name2 == "A. arizonae":
            col2 = "#1469A7"      
        elif args.species_name2 == "A. neomexicanus":
            col2 = "#688e26"
        else:
            col2 = "#b99666"
    
    if args.input_fasta3 is not None:
        if args.plot_color3 is not None:
            col3 = args.plot_color3
        elif args.species_name3 == "A. marmoratus":
            col3 = "#b99666"
        elif args.species_name3 == "A. arizonae":
            col3 = "#1469A7"      
        elif args.species_name3 == "A. neomexicanus":
            col3 = "#688e26"
        else:
            col3 = "#b99666"
    

    label1 = str(args.species_name1) if args.species_name1 else "species1"
    species_order.append(label1)
    jobs.append((label1, args.input_fasta1, args.input_pep1, col1))
    color_map[label1] = col1

    # species 2 (optional – do NOT require a name)
    if args.input_fasta2 and args.input_pep2:
        label2 = str(args.species_name2) if args.species_name2 else "species2"
        species_order.append(label2)
        jobs.append((label2, args.input_fasta2, args.input_pep2, col2))
        color_map[label2] = col2

    # species 3 (optional)
    if args.input_fasta3 and args.input_pep3:
        label3 = str(args.species_name3) if args.species_name3 else "species3"
        species_order.append(label3)
        jobs.append((label3, args.input_fasta3, args.input_pep3, col3))
        color_map[label3] = col3
        
    # summarize each species
    rows = []
    color_map = {}
    for sp, fa, pep, col in jobs:
        rows.append(summarize_species_row(sp, fa, pep))
        color_map[sp] = col

    # table --> one row per species with --> columns = categories
    wide_df = pd.DataFrame(rows)
    # long form table
    long_df = wide_df.melt(
        id_vars=["species"],
        value_vars=["total_transcripts", "with_orf", "complete", "5partial", "3partial", "internal"],
        var_name="category",
        value_name="count"
    )

    # pivot: rows = category, columns = species, values = count
    pivot = long_df.pivot(index="category", columns="species", values="count").fillna(0).astype(int)
    pivot = pivot[species_order]
    
    percent_pivot = pivot.divide(pivot.loc["total_transcripts"], axis=1) * 100
    percent_pivot.loc["total_transcripts"] = 100.0
    percent_pivot = percent_pivot.round(2)

    # save
    tab_1 = os.path.join(args.output_path, "5_orf_statistics_pivot.csv")
    tab_2 = os.path.join(args.output_path, "5_orf_statistics_long.csv")
    pivot.to_csv(tab_1)
    long_df.to_csv(tab_2, index=False)
    print(f"Tables saved to: {tab_1} & {tab_2}")

    # plot grouped bars from pivot
    plot_grouped_from_pivot(pivot, color_map, args.output_path)

    print("Pivot table (rows = category, columns = species, values = count):")
    print(pivot)
    print(percent_pivot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # required species
    parser.add_argument("input_fasta1", help="Path to Transcriptome FASTA file")
    parser.add_argument("input_pep1", help="Path to TransDecoder PEP file")
    parser.add_argument("-species_name1", help="Name of the species")
    # optional user inputs
    parser.add_argument("-plot_color1", help="Color for plotting in hex code")
    # optional species 2
    parser.add_argument("-input_fasta2", help="Path to Transcriptome FASTA file")
    parser.add_argument("-input_pep2", help="Path to TransDecoder PEP file")
    parser.add_argument("-species_name2", help="Name of the species")
    parser.add_argument("-plot_color2", help="Color for plotting in hex code, default") 
    # optional species 3
    parser.add_argument("-input_fasta3", help="Path to Transcriptome FASTA file")
    parser.add_argument("-input_pep3", help="Path to TransDecoder PEP file")
    parser.add_argument("-species_name3", help="Name of the species")
    parser.add_argument("-plot_color3", help="Color for plotting in hex code, default")  
    #optional file prefix
    parser.add_argument("-output_path", help="Output path/folder") 

    args = parser.parse_args()


    main()
