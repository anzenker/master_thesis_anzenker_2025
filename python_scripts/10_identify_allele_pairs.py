#!/usr/bin/env python3
import os, re, tempfile
from time import time
import pandas as pd
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline

# ---------------------------
# LOAD INPUT FILES
# ---------------------------
def load_pep_dict(pep_fasta):
    
    kept_dict = {}
    n_total = n_p1 = 0
    
    # iterate through pep file and kepp only .p1 seq IDs
    for rec in SeqIO.parse(pep_fasta, "fasta"):
        n_total += 1
        rid = rec.id.split()[0]
        if not rid.endswith(".p1"):
            continue
        n_p1 += 1
        seq = str(rec.seq).upper().replace("*", "")
        if seq:
            kept_dict[rid] = seq  
            
    print(f"PEP file: {os.path.basename(pep_fasta)};\n total={n_total}, .p1={n_p1}, kept={len(kept_dict)}\n")
    return kept_dict

def gff_to_transcript_table(gff_path, suffix):

    # GFF file columns
    cols = [f"{suffix}_chr", f"{suffix}_source", f"{suffix}_feature",
            f"{suffix}_start", f"{suffix}_end", f"{suffix}_score",
            f"{suffix}_strand", f"{suffix}_frame", f"{suffix}_attribute"]
    
    df = pd.read_csv(gff_path, sep="\t", header=None, names=cols, comment="#")

    # only keep transcript rows
    transcripts = df[df[f"{suffix}_feature"] == "transcript"].copy()

    # get info from attributes
    transcripts[f"{suffix}_transcript_id"] = transcripts[f"{suffix}_attribute"].str.extract(r'\bID=([^;]+)', expand=False)
    transcripts[f"{suffix}_gene_id"]       = transcripts[f"{suffix}_attribute"].str.extract(r'\b(?:geneID|gene_id)=([^;]+)', expand=False)
    # UniProt annotation --> blast_id
    transcripts[f"{suffix}_blast_id"]      = transcripts[f"{suffix}_attribute"].str.extract(r'\bblast_id=([^;]+)', expand=False) 
    
    # only keep specified cols
    keep = [f"{suffix}_transcript_id", f"{suffix}_gene_id", f"{suffix}_chr", f"{suffix}_blast_id"]
    out = transcripts.loc[:, keep].drop_duplicates()
    
    print(f"GFF file: {os.path.basename(gff_path)}; \n {suffix} transcripts loaded: {len(out)} \n")
    
    return out

def blastm6_to_df(path, suffix1, suffix2):

    s1 = suffix1
    s2 = suffix2

    # col names outfmt6 file
    cols = [
        f"{s1}_transcript_id", f"{s2}_transcript_id",
        f"pident_{s1}2{s2}", f"length_{s1}2{s2}",
        f"mismatch_{s1}2{s2}", f"gapopen_{s1}2{s2}",
        f"qstart_{s1}", f"qend_{s1}", f"sstart_{s2}", f"send_{s2}",
        f"evalue_{s1}2{s2}", f"bitscore_{s1}2{s2}"
    ]

    df = pd.read_csv(path, sep="\t", names=cols)
    
    # only keep specified columns
    keep = [f"{s1}_transcript_id", f"{s2}_transcript_id", f"pident_{s1}2{s2}", f"length_{s1}2{s2}",
            f"evalue_{s1}2{s2}", f"bitscore_{s1}2{s2}"]
    out = df.loc[:, keep]
    
    print(f"BLAST file ({s1}2{s2}): {os.path.basename(path)};\n rows: {len(out)} \n")
    
    return out

# ---------------------------
# IDENTIFY RBH best hit pairs
# ---------------------------
def rbh_from_two(blast1, blast2, suffix1, suffix2):

    s1 = suffix1
    s2 = suffix2

    best_1 = blast1.sort_values(f"bitscore_{s1}2{s2}", ascending=False).drop_duplicates(f"{s1}_transcript_id")
    best_2 = blast2.sort_values(f"bitscore_{s2}2{s1}", ascending=False).drop_duplicates(f"{s2}_transcript_id")
    # remove columns to not merge duplicated rows between blast 1 & blast 2
    best_2 = best_2[[f"{s1}_transcript_id", f"{s2}_transcript_id"]]
    
    rbh = (best_1.merge(best_2, on=[f"{s1}_transcript_id", f"{s2}_transcript_id"], how="inner").drop_duplicates())
    
    print(f"all RBH pairs ({s1}&{s2}): {len(best_1)} best-{s1}, {len(best_2)} best-{s2};\n FINAL RBH pairs merged: {len(rbh)}")
    
    return rbh

# ---------------------------
# EMBOSS needle
# ---------------------------
def run_needle(seq_a_id, seq_b_id, seqA, seqB, out_dir, gapopen: float = 10.0, gapextend: float = 0.5):

    os.makedirs(out_dir, exist_ok=True)
    sa, sb = seq_a_id, seq_b_id

    fa_aln = os.path.join(out_dir, f"{sa}__vs__{sb}.faaln")
    pair_txt = os.path.join(out_dir, f"{sa}__vs__{sb}.pair")

    # Write temporary FASTAs for this pair
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as qa,\
        tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tb:
        qa.write(f">{seq_a_id}\n{seqA}\n")
        tb.write(f">{seq_b_id}\n{seqB}\n")
        qpath, tpath = qa.name, tb.name

    try:
        # FASTA alignment output
        NeedleCommandline(cmd="needle", asequence=qpath, bsequence=tpath,
                        gapopen=gapopen, gapextend=gapextend,
                        aformat="fasta", auto=True, outfile=fa_aln)()

        # output readable alignment .pair FILE
        NeedleCommandline(cmd="needle", asequence=qpath, bsequence=tpath,
                        gapopen=gapopen, gapextend=gapextend,
                        aformat="pair", auto=True, outfile=pair_txt)()

        # go through FASTA alignment file
        records = list(SeqIO.parse(fa_aln, "fasta"))
        seq1, seq2 = str(records[0].seq), str(records[1].seq)
        # exclud gaps "-"
        both = [(x, y) for x, y in zip(seq1, seq2) if x != "-" and y != "-"]
        matches = sum(x == y for x, y in both)
        pid = 100.0 * matches / len(both)
        cov1 = sum(x != "-" for x in seq1) / max(1, len(seqA))
        cov2 = sum(y != "-" for y in seq2) / max(1, len(seqB))

        return {"pid_overlap": round(pid, 2),
                "cov_a": round(cov1, 3), "cov_b": round(cov2, 3),
                "aln_len": len(seq1), "aln_fasta": fa_aln, "aln_pair": pair_txt}
    finally:
        for p in (qpath, tpath):
            try: os.remove(p)
            except: pass

def map_tx_to_prot_id(transcript_id):
    # add .p1 to transcript ID to make matching preotin ID easier
    return f"{transcript_id}.p1"

def protein_align_pairs_with_needle(rbh_df, suffix1, suffix2, pep1_dict, pep2_dict, out_dir, gapopen = 10.0, gapextend = 0.5,):

    s1 = suffix1
    s2 = suffix2

    out_dir = os.path.join(out_dir, "protein_alignments")
    os.makedirs(out_dir, exist_ok=True)

    r = rbh_df.copy()
    r[f"{s1}_prot_id"] = r[f"{s1}_transcript_id"].astype(str).map(map_tx_to_prot_id)
    r[f"{s2}_prot_id"] = r[f"{s2}_transcript_id"].astype(str).map(map_tx_to_prot_id)

    # keep only pairs where both proteins exist
    mask = (r[f"{s1}_prot_id"].isin(pep1_dict) & r[f"{s2}_prot_id"].isin(pep2_dict))
    r = r.loc[mask].reset_index(drop=True)
    print(f"[needle {s1}|{s2}] will align {len(r)} pairs")

    results = []
    aligned = skipped = 0
    t0 = time()
    # align each pair and get alignment stats
    for i, row in r.iterrows():
        a_id = row[f"{s1}_prot_id"]
        b_id = row[f"{s2}_prot_id"]
        stats = run_needle(a_id, b_id, pep1_dict[a_id], pep2_dict[b_id], out_dir=out_dir, gapopen=gapopen, gapextend=gapextend)
        
        aligned += 1
        
        results.append({
                f"{s1}_transcript_id": row[f"{s1}_transcript_id"],
                f"{s2}_transcript_id": row[f"{s2}_transcript_id"],
                f"{s1}_prot_id": a_id, f"{s2}_prot_id": b_id,
                "prot_pid_overlap": stats["pid_overlap"],
                "prot_cov_1": stats["cov_a"],
                "prot_cov_2": stats["cov_b"],
                "aln_len": stats["aln_len"],
                "aln_pair": stats["aln_pair"],
                "aln_fasta": stats["aln_fasta"],
        })


        if (i + 1) % 100 == 0 or (i + 1) == len(r):
            pct = 100.0 * (i + 1) / len(r)
            print(f"  Progress: {i+1}/{len(r)} ({pct:.1f}%) | aligned={aligned}, skipped={skipped} | {(time()-t0)/60:.1f} min")

    print(f"NEEDLE ALIGNMENT {s1}|{s2}] done: aligned={aligned}, skipped={skipped}, total={len(r)}")
    return pd.DataFrame(results)


# ---------------------------
# RUN ALL STEPS 
# ---------------------------
def run_pipeline(name, suffix1, suffix2, paths_dict, synteny_table, syn_cols_for_sides, output_dir):

    s1 = suffix1
    s2 = suffix2

    out_dir = os.path.join(output_dir, name)
    os.makedirs(out_dir, exist_ok=True)
    
    #==================================================================
    # 1) load input files
    anno1 = gff_to_transcript_table(paths_dict["gff1"], s1)
    anno2 = gff_to_transcript_table(paths_dict["gff2"], s2)

    rbh_blast_12 = blastm6_to_df(paths_dict["blast_s1_to_s2"], s1, s2)
    rbh_blast_21 = blastm6_to_df(paths_dict["blast_s2_to_s1"], s2, s1)

    syntenic_df = pd.read_csv(synteny_table, sep="\t", names=[f"{s1}_chr", f"{s1}_chr_no", f"{s2}_chr", f"{s2}_chr_no"])

    pep1 = load_pep_dict(paths_dict["pep1"])
    pep2 = load_pep_dict(paths_dict["pep2"])

    #==================================================================
    # 2) FILTER best hit RBH from both directions (s1 -> s2) & (s2 -> s1)
    rbh = rbh_from_two(rbh_blast_12, rbh_blast_21, s1, s2)
    print(f"[{name}] after RBH: {len(rbh)}")
    rbh.to_csv(os.path.join(out_dir, f"{name}_rbh.tsv"), sep="\t", index=False)

    #==================================================================
    # 3) FILTER for shared UniProt filter
    rbh_with_anno_one = rbh.merge(anno1, on=f"{s1}_transcript_id", how="left")
    rbh_with_anno = rbh_with_anno_one.merge(anno2, on=f"{s2}_transcript_id", how="left")

    rbh_with_anno["shared_blast_id"] = (rbh_with_anno[f"{s1}_blast_id"].notna() 
                                        & rbh_with_anno[f"{s2}_blast_id"].notna()
                                        & (rbh_with_anno[f"{s1}_blast_id"] == rbh_with_anno[f"{s2}_blast_id"]))

    rbh_shared = rbh_with_anno[rbh_with_anno["shared_blast_id"]].copy()
    print(f"FILTER shared UniProt - kept RBH pairs: {len(rbh_shared)}")
    rbh_shared.to_csv(os.path.join(out_dir, f"{name}_rbh_shared_uniprot.tsv"), sep="\t", index=False)

    #==================================================================
    # 4) sanity check - FILTER for syntenic chromsome origin
    chr_1 = f"{s1}_chr"
    chr_2 = f"{s2}_chr"
    left  = rbh_shared.copy()
    right = syntenic_df[[chr_1, chr_2]].copy()
    # only keep matches between dfs
    rbh_synt = left.merge(right, on=[chr_1, chr_2], how="inner")

    print(f"{name} after synteny: {len(rbh_synt)}")
    rbh_synt.to_csv(os.path.join(out_dir, f"{name}_rbh_shared_uniprot_syntenic.tsv"), sep="\t", index=False)

    #==================================================================
    # 5) Global alignment --> protein alignments (needle) using .p1 mapping
    print(f"[{s1}] PEP entries: {len(pep1)} | [{s2}] PEP entries: {len(pep2)}")

    aln_df = protein_align_pairs_with_needle(rbh_synt, s1, s2, pep1, pep2, out_dir)
    aln_df.to_csv(os.path.join(out_dir, f"{name}_protein_alignment_summary.tsv"), sep="\t", index=False)

    #==================================================================
    # 6) output check - counts in between filtering steps
    print(f"{name} FINAL:")
    print(f"  RBH pairs total                 : {len(rbh)}")
    print(f"  Shared UniProt (both annotated) : {len(rbh_shared)}")
    print(f"  Syntenic                         : {len(rbh_synt)}")
    print(f"  Protein-alignable (.p1 present) : {len(aln_df)}")

    return {"rbh": rbh, "rbh_shared": rbh_shared, "rbh_syntenic": rbh_synt, "protein_alignments": aln_df, "out_dir": out_dir}

# ---------------------------
# Example usage
# ---------------------------
if __name__ == "__main__":
    # INPUT synteny table
    # cols: Ari_Chr_X   Marm_Chr_X
    synteny_tsv = "synteny_table.txt"
    
    # PATH DICT --> comparing A. arizonae to A. marmoratus
    id_allele_pairs_ari_marm = {
        "gff1": "/RESULTS_Aari/9_GFF_annotation_files/stringtie2_ASP_ARI_transcriptome_isoforms_inherit_anno.gff",
        "gff2": "/RESULTS_Amarm/9_GFF_annotation_files/stringtie2_ASP_MARM_transcriptome_isoforms_inherit_anno.gff",
        "blast_s1_to_s2": "blastn_ari2marmDB_transcriptome.outfmt6",
        "blast_s2_to_s1": "blastn_marm2ariDB_transcriptome.outfmt6",
        "pep1": "/RESULTS_Aari/5_frame_selection/stringtie2_ASP_ARI_transcriptome_renamed.fasta.transdecoder_dir/longest_orfs.pep",
        "pep2": "/RESULTS_Amarm/5_frame_selection/stringtie2_ASP_MARM_transcriptome_renamed.fasta.transdecoder_dir//longest_orfs.pep"
    }

    # PATH DICT --> A. neomexicanus --> comparing A. arizonae transcripts to A. marmoratus transcripts
    id_allele_pairs_neo = {
        "gff1": "/RESULTS_Aneo/9_GFF_annotation_files/stringtie2_ASP_NEO_transcriptome_isoforms_inherit_anno.gff",  # neo_ari GFF
        "gff2": "/RESULTS_Aneo/9_GFF_annotation_files/stringtie2_ASP_NEO_transcriptome_isoforms_inherit_anno.gff",   # neo_marm GFF (same file if combined, OK)
        "blast_s1_to_s2": "/RESULTS_Aneo/blastn_neo_ari2marmDB_transcriptome.outfmt6",
        "blast_s2_to_s1": "/RESULTS_Aneo/blastn_neo_marm2ariDB_transcriptome.outfmt6",
        "pep1": "/RESULTS_Aneo/5_frame_selection/stringtie2_ASP_NEO_transcriptome_renamed.fasta.transdecoder_dir/longest_orfs.pep",
        "pep2": "/RESULTS_Aneo/5_frame_selection/stringtie2_ASP_NEO_transcriptome_renamed.fasta.transdecoder_dir/longest_orfs.pep"
        }

    out_path = "10_identifying_allele_pairs_output"
    os.makedirs(out_path, exist_ok=True)

    # ARI against MARM
    run_pipeline("ari_vs_marm", "ari", "marm", id_allele_pairs_ari_marm, synteny_tsv, ("Ari_Chr","Marm_Chr"), out_path)

    # NEO_ARI against NEO_MARM
    run_pipeline("neo_ari_vs_neo_marm", "neo_ari", "neo_marm", id_allele_pairs_neo, synteny_tsv, ("Ari_Chr","Marm_Chr"), out_path)