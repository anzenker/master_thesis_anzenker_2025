#!/usr/bin/env python3
import argparse
import re

def clean_id(qid: str) -> str:
    """Drop trailing .p# on peptide IDs (e.g., foo.p1 -> foo)."""
    return re.sub(r"\.p\d+$", "", qid.strip())

def load_eggnog_map(path: str) -> dict:
    """
    Map query -> {'eggnog_taxa': 'xxxx', 'eggnog_id': 'YYY'}
    Accepts .annotations or .hits. Uses first two tab columns.
    """
    m = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            q, seed = line.rstrip("\n").split("\t", 2)[:2]
            q = clean_id(q)
            taxa, eid = (seed.split(".", 1) + [""])[:2] if "." in seed else ("", seed)
            m[q] = {"eggnog_taxa": taxa, "eggnog_id": eid}
    return m

def load_blast_map(path: str) -> dict:
    """
    Map query -> {'blast_db': 'sp', 'blast_id': 'Q9XXX', 'blast_name': 'NAME'}
    Assumes outfmt6 col1 looks like db|ACC|NAME (fallbacks if not).
    """
    m = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            q, s = line.rstrip("\n").split("\t", 2)[:2]
            q = clean_id(q)
            parts = s.split("|")
            if len(parts) >= 3:
                db, acc, name = parts[0], parts[1], parts[2]
            else:
                db, acc, name = "", s, ""
            m[q] = {"blast_db": db, "blast_id": acc, "blast_name": name}
    return m

# ----------------- GFF3 attr parsing -----------------

def parse_gff_attrs(attr_str: str) -> dict:
    """
    Parse GFF3 attributes 'key=value;key=value;...' into dict (no URL-decoding).
    """
    attrs = {}
    for field in attr_str.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            attrs[k] = v
        else:
            attrs[field] = ""
    return attrs

def gff_attrs_to_str(attrs: dict) -> str:
    # deterministic order
    return ";".join([f"{k}={v}" if v != "" else k for k, v in sorted(attrs.items())])


def annotate_gff(in_gff: str, out_gff: str, eggnog_map: dict | None, blast_map: dict | None, annotate_genes_too: bool = False):
    """
    Read GFF3, add EggNOG / BLAST fields to selected features, write GFF3.
    - Default: annotate only feature == 'transcript' (using 'ID')
    - If annotate_genes_too: also annotate feature == 'gene' (using 'ID')
    If your transcripts are identified via 'transcript_id' instead of 'ID',
    this script will also try 'transcript_id' as a fallback.
    """
    eggnog_map = eggnog_map or {}
    blast_map  = blast_map  or {}

    with open(in_gff) as fin, open(out_gff, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line); continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line); continue

            feature = cols[2]
            attrs = parse_gff_attrs(cols[8])

            # Which rows to annotate?
            if feature == "transcript":
                key_id = attrs.get("ID") or attrs.get("transcript_id")
            elif feature == "gene" and annotate_genes_too:
                key_id = attrs.get("ID")
            else:
                fout.write(line); continue

            if not key_id:
                fout.write(line); continue

            base_id = clean_id(key_id)

            # attach eggNOG (do not overwrite existing)
            if base_id in eggnog_map:
                e = eggnog_map[base_id]
                attrs.setdefault("eggnog_taxa", e["eggnog_taxa"])
                attrs.setdefault("eggnog_id",   e["eggnog_id"])

            # attach BLAST (do not overwrite existing)
            if base_id in blast_map:
                b = blast_map[base_id]
                attrs.setdefault("blast_db",   b["blast_db"])
                attrs.setdefault("blast_id",   b["blast_id"])
                attrs.setdefault("blast_name", b["blast_name"])

            cols[8] = gff_attrs_to_str(attrs)
            fout.write("\t".join(cols) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Append eggNOG and/or UniProt BLAST annotations to a GFF3."
    )
    parser.add_argument("--gff",  required=True, help="Input GTF")
    parser.add_argument("--out",  required=True, help="Output GTF")
    parser.add_argument("--eggnog", help="eggNOG .annotations/.hits (optional)")
    parser.add_argument("--blast",  help="BLAST outfmt6 vs UniProt (optional)")
    parser.add_argument("--features", help="One feature or comma-separated feature types to annotate. Example: transcript,gene")
    args = parser.parse_args()

    if not args.eggnog and not args.blast:
        parser.error("Provide at least one of --eggnog or --blast.")

    features = [f.strip() for f in args.features.split(",") if f.strip()]
    if not features:
        parser.error("Empty --features. Example: --features transcript,gene")
    
    if "gene" in features:
        annotate_genes_too = True
    else: 
        annotate_genes_too = False

    eggnog_map = load_eggnog_map(args.eggnog) if args.eggnog else None
    blast_map  = load_blast_map(args.blast)   if args.blast  else None

    annotate_gff(args.gff, args.out, eggnog_map, blast_map, annotate_genes_too)

if __name__ == "__main__":
    main()