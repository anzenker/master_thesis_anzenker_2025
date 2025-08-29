#!/usr/bin/env python3

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Pick canonical transcript IDs with max coverage per gene.")
    parser.add_argument("input_tsv", help="Input TSV file with columns: gene_id, transcript_id, coverage")
    parser.add_argument("output_txt", help="Output TXT file for canonical transcript IDs")
    args = parser.parse_args()

    # read input file
    df = pd.read_csv(args.input_tsv, sep='\t', header=None, names=["gene_id", "transcript_id", "coverage"])

    # make sure coverage is float
    df["coverage"] = df["coverage"].astype(float)

    # keep isoform per gene with max coverage
    canonical_df = df.loc[df.groupby("gene_id")["coverage"].idxmax()]

    # save transcript IDs to output file
    canonical_df["transcript_id"].to_csv(args.output_txt, sep="\t", header=False, index=False)

if __name__ == "__main__":
    main()