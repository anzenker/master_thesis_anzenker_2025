#!/usr/bin/env python3
import argparse, textwrap, os, gzip
import pandas as pd
import matplotlib.pyplot as plt

# ---------- helpers ----------

SPECIES_COLORS = {
    "A. marmoratus":  "#b99666",
    "A. arizonae":    "#1469A7",
    "A. neomexicanus":"#688e26",
}

def count_fasta_headers(path: str) -> int:
    #count sequences by '>' header
    open_fa = gzip.open if path.endswith(".gz") else open
    with open_fa(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith(">"))

def summarize_species_counts(species: str, fasta: str, canonical_fasta: str) -> pd.DataFrame:
    # df count total and canonical
    total_transcripts = count_fasta_headers(fasta)
    total_canonical   = count_fasta_headers(canonical_fasta)

    df = pd.DataFrame(
        [
            ["Total Count Transcripts", total_transcripts, species],
            ["Total Count Canonical",   total_canonical,   species],
        ],
        columns=["Total", "Count", "Species"]
    )
    return df

def default_color_for(species: str) -> str:
    return SPECIES_COLORS.get(species.strip(), "#C79FEF")

def plot_transcript_counts(df: pd.DataFrame, output_path: str, color_map: dict):
    # order species as they appear
    species_order = list(dict.fromkeys(df["Species"].tolist()))

    df["Total"] = pd.Categorical(
        df["Total"],
        categories=["Total Count Transcripts", "Total Count Canonical"],
        ordered=True
    )
    df["Species"] = pd.Categorical(df["Species"], categories=species_order, ordered=True)

    pivot_df = df.pivot(index="Total", columns="Species", values="Count").fillna(0).astype(int)
    colors = [color_map.get(sp, default_color_for(sp)) for sp in pivot_df.columns]

    ax = pivot_df.plot(kind="bar", width=0.5, edgecolor="white", color=colors, figsize=(9, 6))
    ax.set_ylabel("Count of Transcripts")
    ax.set_xlabel("")
    ax.set_title("Total vs Canonical Transcript Counts per Species")
    plt.xticks(rotation=0)
    ax.legend(title="Species", frameon=False)

    species = list(pivot_df.columns)          # columns = species
    cats     = list(pivot_df.index)           # rows = categories
    totals_by_species = pivot_df.loc["Total Count Transcripts"].to_dict()  # {species: total}

    
    # annotate bars --> outcomment if not wanted
    ymax = pivot_df.max().max() * 1.10
    ax.set_ylim(0, ymax)
    for col_idx, container in enumerate(ax.containers):
        sp = species[col_idx]
        denom = float(totals_by_species.get(sp, 0)) or 1.0  # avoid division by 0

        for j, bar in enumerate(container):                 # j = category index
            count = bar.get_height()
            # Percent of that species' total transcripts
            pct = 100.0 * count / denom if denom else 0.0

            ax.text(
                bar.get_x() + bar.get_width()/2,
                count + ymax*0.01,
                f"{pct:.1f}%\n{int(count)}",
                ha="center", va="bottom", fontsize=9
        )
    plt.tight_layout()

    # save
    fig = os.path.join(output_path, "4_plot_total_vs_canonical_counts.png")
    plt.savefig(fig, dpi=500)
    print(f"Plot saved to: {fig}")

    plt.close()


def main():
    p = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            ...
        """),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # species 1 (required FASTA, optional name)
    p.add_argument("fasta1", help="Transcriptome FASTA (species 1)")
    p.add_argument("canonical1", help="Canonical transcript FASTA (species 1)")
    p.add_argument("--species1", default=None, help="Species name 1 (optional)")
    p.add_argument("--color1", default=None, help="Hex color for species 1 (optional)")

    # species 2 (optional)
    p.add_argument("--fasta2", default=None, help="Transcriptome FASTA (species 2)")
    p.add_argument("--canonical2", default=None, help="Canonical transcript FASTA (species 2)")
    p.add_argument("--species2", default=None, help="Species name 2 (optional)")
    p.add_argument("--color2", default=None, help="Hex color for species 2 (optional)")

    # species 3 (optional)
    p.add_argument("--fasta3", default=None, help="Transcriptome FASTA (species 3)")
    p.add_argument("--canonical3", default=None, help="Canonical transcript FASTA (species 3)")
    p.add_argument("--species3", default=None, help="Species name 3 (optional)")
    p.add_argument("--color3", default=None, help="Hex color for species 3 (optional)")

    p.add_argument("--output_path", help="Output path/folder")
    args = p.parse_args()

    data, color_map = [], {}

    # s1
    sp1 = args.species1 or "Species1"
    data.append((sp1, args.fasta1, args.canonical1))
    color_map[sp1] = args.color1 or default_color_for(sp1)

    # s2
    if args.fasta2 and args.canonical2:
        sp2 = args.species2 or "Species2"
        data.append((sp2, args.fasta2, args.canonical2))
        color_map[sp2] = args.color2 or default_color_for(sp2)

    # s3
    if args.fasta3 and args.canonical3:
        sp3 = args.species3 or "Species3"
        data.append((sp3, args.fasta3, args.canonical3))
        color_map[sp3] = args.color3 or default_color_for(sp3)

    # put data into df
    species = [summarize_species_counts(species, fasta_all, canonical) for species, fasta_all, canonical in data]
    df = pd.concat(species, ignore_index=True)

    # save
    path_tab_1 = os.path.join(args.output_path, "4_total_vs_canonical_counts.csv")
    df.to_csv(path_tab_1, index=False)
    print(f"Saved table: {tab_1}")

    plot_transcript_counts(df, args.output_path, color_map)

if __name__ == "__main__":
    main()
