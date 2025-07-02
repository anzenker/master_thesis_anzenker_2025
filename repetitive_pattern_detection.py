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

def find_repetitive_patterns(sequence):
    patterns = {
        "repetitive_AC": r"(AC|CA|ACA|AAC|ACC|CAC){10,}",
        "repetitive_GT": r"(GT|TG|TGT|TGG|GTT|TGT){10,}",
        "repetitive_AG": r"(GA|AG|AAGA|AGA|GAA|GAG){10,}",
        "repetitive_AT": r"(AT|TA|AAT|AAAT|TAT|ATA|TAA){10,}",
        "repetitive_CG": r"(GC|CG|CGC|CGG|GCC|CGC){10,}",
        "repetitive_CT": r"(TC|CT|CTC|TCC|CTT|CTC){10,}",
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

        log_print(f"Starting repetitive pattern detection for {input_fasta}", log)

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
        hit_read_ids = set()
        total_reads = 0

        log_print("Reading input FASTA / FASTQ file...", log)
        with open_func(input_fasta, "rt") as handle:
            for record in SeqIO.parse(handle, seq_format):
                total_reads += 1
                read_id = record.id
                sequence = str(record.seq)
                matches = find_repetitive_patterns(sequence)
                for pattern_type, start, end, length in matches:
                    pattern_hits.append([read_id, len(sequence), pattern_type, start, end, length])
                    hit_read_ids.add(read_id)

        log_print(f"Scanned {total_reads} reads. {len(hit_read_ids)} reads contained repetitive patterns.", log)
        log_print(f"{(len(hit_read_ids) / total_reads) * 100} % of reads contain a repetitive pattern.", log)

        # Save read IDs to txt
        log_print(f"Writing read IDs with patterns to: {output_ids_path}", log)
        with open(output_ids_path, 'w') as f:
            for read_id in sorted(hit_read_ids):
                f.write(f"{read_id}\n")

        # Save detailed TSV
        log_print(f"Writing detailed pattern info to: {output_tsv_path}", log)
        df = pd.DataFrame(pattern_hits, columns=["read_id", "seq_length", "pattern_type", "start", "end", "pattern_length"])
        df.to_csv(output_tsv_path, sep='\t', index=False)

        elapsed_time = time.time() - start_time
        log_print(f"Finished in {elapsed_time:.2f} seconds.", log)
        log_print(f"Log saved to: {log_path}", log)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect repetitive patterns (>10 nt) in a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA / FASTA.gz / FASTQ / FASTQ.gz file")
    parser.add_argument("output_ids", help="Output TXT file with read IDs containing repetitive patterns")
    parser.add_argument("output_tsv", help="Output TSV file with pattern info per read")
    args = parser.parse_args()

    # Automatically name the log file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{os.path.splitext(args.output_tsv)[0]}_log_{timestamp}.txt"

    process_fasta(args.input_fasta, args.output_ids, args.output_tsv, log_file)