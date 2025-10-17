#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from datetime import datetime
import gzip
import os
import pandas as pd
import re
import sys
import time

def open_func(input_file):
    #open gzip file
    return gzip.open if input_file.endswith(".gz") else open

def log_print(text, log_handle):
    #logging file to keep track of run
    print(text)

    log_handle.write(text + "\n")

def find_repetitive_patterns(sequence):
    matches = []

    # defined patterns to detect, ANYWHERE in the sequence
    patterns = {
        "repetitive_AC": r"(AC|CA|ACA|AAC|ACC|CAC){10,}",
        "repetitive_GT": r"(GT|TG|TGT|TGG|GTT|TGT){10,}",
        "repetitive_AG": r"(GA|AG|AAGA|AGA|GAA|GAG){10,}",
        "repetitive_AT": r"(AT|TA|AAT|AAAT|TAT|ATA|TAA){10,}",
        "repetitive_CG": r"(GC|CG|CGC|CGG|GCC|CGC){10,}",
        "repetitive_CT": r"(TC|CT|CTC|TCC|CTT|CTC){10,}",
    }

    # iterates over each defined pattern and 
    # checks whether it is found in the seqeunce
    for pattern_name, pattern in patterns.items():
        # iterates with an index
        for match in re.finditer(pattern, sequence):
            start = match.start()
            end = match.end()
            length = end - start
            # only detect & document patterns with a length > 10
            if length > 10:
                matches.append((pattern_name, start, end, length))

    return matches

def process_sequences(input_file, output_ids_path, output_tsv_path, log_path):

    pattern_hits = []
    hit_read_ids = set()
    total_reads = 0
    
    start_time = time.time()

    # open log file to document script process --> 'w' = write
    with open(log_path, 'w') as log:

        log_print(f"Starting repetitive pattern detection for {input_file}", log)

        # determine sequence format: fasta or fastq
        if input_file.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
            seq_format = "fasta"
        elif input_file.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            seq_format = "fastq"
        else:
            log_print("Error: Not a valid file format. Please input a .fasta/.fa/.fastq/.fq with optional .gz file.", log)
            # exit script
            sys.exit(1)
        
        log_print(f"Detected input file format: {seq_format}", log)

        log_print("Reading input FASTA / FASTQ file...", log)
        # open input file, 'rt' = read & write
        with open_func(input_file)(input_file, "rt") as handle:
            # handle FASTA/FASTQ file
            for record in SeqIO.parse(handle, seq_format):
                total_reads += 1
                read_id = record.id
                sequence = str(record.seq)
                matches = find_repetitive_patterns(sequence)
                for pattern_type, start, end, length in matches:
                    # save all detected patterns with additional info
                    pattern_hits.append([read_id, len(sequence), pattern_type, start, end, length])
                    # save read ids of seqeunces with a detected pattern
                    hit_read_ids.add(read_id)

        log_print(f"Scanned {total_reads} reads. {len(hit_read_ids)} reads contained repetitive patterns.", log)
        log_print(f"{(len(hit_read_ids) / total_reads) * 100} % of reads contain a repetitive pattern.", log)

        # save read ids of seqeunces with a detected pattern to a TXT file
        log_print(f"Writing read IDs with patterns to: {output_ids_path}", log)
        with open(output_ids_path, 'w') as txt_file:
            for read_id in sorted(hit_read_ids):
                txt_file.write(f"{read_id}\n")

        # save detailed detected pattern information in a TSV file
        log_print(f"Writing detailed pattern info to: {output_tsv_path}", log)
        pattern_df = pd.DataFrame(pattern_hits, columns=["read_id", "seq_length", "pattern_type", "start", "end", "pattern_length"])
        pattern_df.to_csv(output_tsv_path, sep='\t', index=False)

        # determine & save script run time
        elapsed_time = time.time() - start_time
        log_print(f"Finished in {elapsed_time:.2f} seconds.", log)
        log_print(f"Log saved to: {log_path}", log)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Detect repetitive patterns (>10 nt) in a FASTA/FASTQ file.")
    parser.add_argument("input_file", help="Input FASTA / FASTA.gz / FASTQ / FASTQ.gz file")
    args = parser.parse_args()

    # input file basename
    basename_file = os.path.basename(args.input_file)

    # automatically name the log file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{basename_file}_log_repetitive_pattern_detection_{timestamp}.txt"

    # automatically name TXT and TSV output file according to input_file
    output_read_ids_txt = f"{basename_file}_read_ids_with repetitive_pattern.txt"
    output_pattern_info_tsv = f"{basename_file}_repetitive_pattern_info.tsv"

    # run script with user input
    process_sequences(args.input_file, output_read_ids_txt, output_pattern_info_tsv, log_file)
