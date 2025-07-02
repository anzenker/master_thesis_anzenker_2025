#!/usr/bin/env python3

import re
import gzip
import argparse
import os
import sys
import time
from Bio import SeqIO
import pandas as pd
from datetime import datetime

def find_polya_patterns(sequence):
    patterns = {
        "polyA_start": r"^A{10,}",
        "polyA_middle": r"A{15,}",
        "polyA_end": r"A{10,}$"
    }
    
    matches = []

    for pattern_name, pattern in patterns.items():
        for match in re.finditer(pattern, sequence):
            start = match.start()
            end = match.end()
            length = end - start
            if length > 10:
                matches.append((pattern_name, start, end, length))

    return matches

def log_print(message, log_handle):
    print(message)
    log_handle.write(message + "\n")

def process_fasta(input_fasta, output_ids_path, output_tsv_path, log_path):
    start_time = time.time()
    with open(log_path, 'w') as log:
        log_print(f"Starting polyA pattern detection for {input_fasta}", log)

        open_func = gzip.open if input_fasta.endswith(".gz") else open
        
        # Determine sequence format: fasta or fastq
        if input_fasta.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
            seq_format = "fasta"
        elif input_fasta.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            seq_format = "fastq"
        else:
            log_print("Error: Could not determine file format. Use .fasta/.fa/.fastq/.fq with optional .gz", log)
            sys.exit(1)
        
        pattern_hits = []
        total_reads = 0

        log_print("Reading input FASTA/FASTQ file...", log)

        with open_func(input_fasta, "rt") as handle:
            for record in SeqIO.parse(handle, seq_format):
                total_reads += 1
                read_id = record.id
                sequence = str(record.seq)
                matches = find_polya_patterns(sequence)
                for pattern_type, start, end, length in matches:
                    pattern_hits.append([read_id, len(sequence), pattern_type, start, end, length])

        df = pd.DataFrame(pattern_hits, columns=["read_id", "seq_length", "pattern_type", "start", "end", "pattern_length"])

        # Remove overlapping polyA_middle
        start_end_ids = set(df[df["pattern_type"].isin(["polyA_start", "polyA_end"])]["read_id"])
        df_filtered = df[~((df["pattern_type"] == "polyA_middle") & (df["read_id"].isin(start_end_ids)))]

        # Get read IDs
        unique_reads_with_polya = set(df_filtered["read_id"])
        percent = (len(unique_reads_with_polya) / total_reads) * 100

        log_print(f"Scanned {total_reads} reads.", log)
        log_print(f"Reads with polyA pattern (after filtering): {len(unique_reads_with_polya)} ({percent:.2f}%)", log)

        # Write read IDs
        with open(output_ids_path, 'w') as f:
            for rid in sorted(unique_reads_with_polya):
                f.write(rid + "\n")

        # Write TSV
        df_filtered.to_csv(output_tsv_path, sep='\t', index=False)
        log_print(f"Pattern info written to: {output_tsv_path}", log)
        log_print(f"Read ID list written to: {output_ids_path}", log)

        elapsed = time.time() - start_time
        log_print(f"Finished in {elapsed:.2f} seconds.", log)

def main():
    parser = argparse.ArgumentParser(description="Detect polyA start/middle/end patterns in a FASTA file.")
    parser.add_argument("input_fasta_OR_fastq", help="Input FASTA / FASTA.gz / FASTQ / FASTQ.gz file")
    parser.add_argument("output_ids", help="Output TXT file with read IDs containing polyA patterns")
    parser.add_argument("output_tsv", help="Output TSV file with pattern info per read")
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{os.path.splitext(args.output_tsv)[0]}_log_{timestamp}.txt"

    process_fasta(args.input_fasta, args.output_ids, args.output_tsv, log_file)

if __name__ == "__main__":
    main()