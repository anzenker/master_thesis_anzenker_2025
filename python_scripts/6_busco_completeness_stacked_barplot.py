import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import textwrap
import numpy as np

def busco_row_from_fulltable(busco_df, species_name):
    #mapping BUSCO status to category name
    status_to_category = {
        "Complete":   "Single-Copy Complete",
        "Duplicated": "Duplicated Complete",
        "Fragmented": "Fragmented",
        "Missing":    "Missing"
    }
    # apply mapping to a new column
    busco_df["category"] = busco_df["Status"].map(status_to_category)

    # Count categories in the desired order
    categories_order = ["Single-Copy Complete", "Duplicated Complete", "Fragmented", "Missing"]

    busco_category_counts = busco_df["category"].value_counts().reindex(categories_order, fill_value=0)
    busco_total_count = busco_category_counts.sum()
    busco_category_percents = round((busco_category_counts / busco_total_count) * 100, 2)

    busco_row = pd.DataFrame({"species" : species_name, "category" : categories_order, "percent" : busco_category_percents})

    return busco_row

def plot_busco_stacked_compare(busco_rows, output_path):

    all_busco_df = pd.concat(busco_rows, ignore_index=True)

    categories_order = ["Single-Copy Complete", "Duplicated Complete", "Fragmented", "Missing"]

    # keep category order consistent
    all_busco_df["category"] = pd.Categorical(all_busco_df["category"], categories=categories_order, ordered=True)

    # pivot: rows = species, columns = category, values = percent
    pivot_busco = all_busco_df.pivot(index="species", columns="category", values="percent").reindex(columns=categories_order)

    colors = {
        "Single-Copy Complete": "#4daf4a",
        "Duplicated Complete":  "#a1d99b",
        "Fragmented":           "#ff7f00",
        "Missing":              "#e41a1c",
        }

    ax = pivot_busco.plot(kind="barh", stacked=True, color=[colors[c] for c in pivot_busco.columns], figsize=(11,3), width=0.2)
    ax.set_title("BUSCO Completeness - Vertebrata Dataset (3357 Orthologs)")
    ax.set_xlabel("Percentage of BUSCOs")
    ax.set_ylabel("")
    ax.set_xlim(0, 110)
    ax.legend(title="BUSCO category", bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)


    # Loop over each category's bars
    for i, container in enumerate(ax.containers):
        for bar in container:
            width = bar.get_width()
            x = bar.get_x() + (width / 2)
            y = bar.get_y() + bar.get_height() * (0.55 if i % 2 == 0 else 0.35)
            ax.text(x, y, f"{width:.1f}%", ha="center", va="center", clip_on=True, fontweight='bold')

    plt.tight_layout()
    
    # Save
    fig = os.path.join(output_path, "6_busco_completeness_stacked_barplot.png")
    plt.savefig(fig)
    print(f"Plot saved to: {fig}")
    
    plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=textwrap.dedent('''\
        ... '''),
        usage='use "python %(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    # input file
    parser.add_argument("full_table_tsv_1", help="Path to busco full_table.tsv")
    parser.add_argument("species_name_1", help="Name of species 1")
    # optional species 2
    parser.add_argument("--full_table_tsv_2", help="Path to busco full_table.tsv")
    parser.add_argument("--species_name_1", help="Name of species 2")
    # optional species 3
    parser.add_argument("--full_table_tsv_3", help="Path to busco full_table.tsv")
    parser.add_argument("--species_name_3", help="Name of species 3")

    parser.add_argument("output_path", help="Path to save output plot image to.")
    args = parser.parse_args()

    col_names = ["BUSCO_id", "Status", "Sequence", "Score", "Length", "OrthoDB", "url", "Description"]

    busco_full_table_df_1 = pd.read_csv(args.full_table_tsv_1, sep="\t", comment="#", names=col_names)

    busco_rows = []

    # species 1 (required)
    df1 = pd.read_csv(args.full_table_tsv_1, sep="\t", comment="#", names=col_names)
    label1 = fr"$\it{{{args.species_name_1}}}$"
    busco_rows.append(busco_row_from_fulltable(df1, label1))

    # species 2 (optional; require both file and name)
    if args.full_table_tsv_2 and args.species_name_2:
        df2 = pd.read_csv(args.full_table_tsv_2, sep="\t", comment="#", names=col_names)
        label2 = fr"$\it{{{args.species_name_2}}}$"
        busco_rows.append(busco_row_from_fulltable(df2, label2))

    # species 3 (optional; require both file and name)
    if args.full_table_tsv_3 and args.species_name_3:
        df3 = pd.read_csv(args.full_table_tsv_3, sep="\t", comment="#", names=col_names)
        label3 = fr"$\it{{{args.species_name_3}}}$"
        busco_rows.append(busco_row_from_fulltable(df3, label3))


    plot_busco_stacked_compare(busco_rows, args.output_path)
