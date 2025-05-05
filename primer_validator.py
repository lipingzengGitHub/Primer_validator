#!/usr/bin/env python3

"""
WGS Primer Validation Tool (Amplicon-based WGS)

Usage:
python primer_validator.py --primer_file primers.csv --reference pathogen_genome.fasta --max_mismatch 2
"""

import argparse
import sys
import os
import re
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
from Bio import SeqIO

def calculate_primer_properties(seq):
    gc_content = gc_fraction(seq) * 100
    tm = mt.Tm_NN(seq)
    return len(seq), gc_content, tm

def reverse_complement(seq):
    complement = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(complement)[::-1]

def fuzzy_search(seq, primer, max_mismatch):
    pattern = ''.join([f"[{base.upper().replace('N', 'ATCG')}]" for base in primer])
    matches = re.finditer(f'(?=({pattern}))', seq)
    for match in matches:
        start = match.start()
        target_seq = seq[start:start+len(primer)]
        mismatches = sum(1 for a, b in zip(target_seq, primer) if a != b and b != 'N')
        if mismatches <= max_mismatch:
            return start
    return -1

def in_silico_pcr(forward_primer, reverse_primer, template, max_mismatch):
    fwd_index = fuzzy_search(template, forward_primer, max_mismatch)
    rev_index = fuzzy_search(template, reverse_complement(reverse_primer), max_mismatch)
    if fwd_index != -1 and rev_index != -1 and fwd_index < rev_index:
        product_size = rev_index + len(reverse_primer) - fwd_index
        return True, product_size
    else:
        return False, None

def check_and_build_blast_db(reference_fasta):
    db_files = [reference_fasta + ext for ext in [".nhr", ".nin", ".nsq"]]
    if all(os.path.exists(f) for f in db_files):
        return reference_fasta
    else:
        makeblastdb_cline = NcbimakeblastdbCommandline(input_file=reference_fasta, dbtype="nucl")
        makeblastdb_cline()
        return reference_fasta

def run_blast(primer_seq, blast_db, out_file):
    blastn_cline = NcbiblastnCommandline(query=primer_seq, db=blast_db, evalue=0.01, outfmt=6, out=out_file)
    stdout, stderr = blastn_cline()
    return os.path.getsize(out_file) > 0

def read_primers(primer_file):
    primers = []
    with open(primer_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    primers.append({"Name": parts[0], "F": parts[1], "R": parts[2]})
    return primers

def read_template(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return str(records[0].seq)

def calculate_score(fwd_tm, rev_tm, off_target, amplifiable):
    score = 100
    tm_diff = abs(fwd_tm - rev_tm)
    score -= tm_diff * 2
    if off_target:
        score -= 30
    if amplifiable:
        score += 10
    return max(score, 0)

def main():
    parser = argparse.ArgumentParser(description="RSV WGS primer validation tool (amplicon-based sequencing)")
    parser.add_argument("--primer_file", required=True, help="CSV file with primer pairs (PrimerName,F,R)")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA file")
    parser.add_argument("--max_mismatch", type=int, default=2, help="Maximum mismatches allowed in in silico PCR (default: 2)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    primers = read_primers(args.primer_file)
    template = read_template(args.reference)
    blast_db = check_and_build_blast_db(args.reference)

    results = []

    for p in primers:
        fwd_len, fwd_gc, fwd_tm = calculate_primer_properties(p["F"])
        rev_len, rev_gc, rev_tm = calculate_primer_properties(p["R"])
        amplifiable, product_size = in_silico_pcr(p["F"], p["R"], template, max_mismatch=args.max_mismatch)

        with open("temp_primer.fasta", "w") as f:
            f.write(f">fwd\n{p['F']}\n>rev\n{p['R']}\n")

        off_target = run_blast("temp_primer.fasta", blast_db, "blast_output.txt")
        score = calculate_score(fwd_tm, rev_tm, off_target, amplifiable)

        results.append({
            "Primer Name": p["Name"],
            "Forward Primer": p["F"],
            "Reverse Primer": p["R"],
            "Fwd_len": fwd_len,
            "Fwd_GC%": fwd_gc,
            "Fwd_Tm": fwd_tm,
            "Rev_len": rev_len,
            "Rev_GC%": rev_gc,
            "Rev_Tm": rev_tm,
            "Amplifiable": amplifiable,
            "Product Size": product_size if amplifiable else "N/A",
            "Off Target?": "Yes" if off_target else "No",
            "Score": score
        })

    df = pd.DataFrame(results)
    df = df.sort_values(by="Score", ascending=False)
    df.to_excel("WGS_primer_validation_full.xlsx", index=False)
    print("Primer validation completed! Results saved to WGS_primer_validation_full.xlsx")

if __name__ == "__main__":
    main()


