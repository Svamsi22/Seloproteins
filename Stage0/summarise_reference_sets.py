# summarise_reference_sets.py
# Prints final summary of both sets before DB building

import argparse
from Bio import SeqIO
import os
from collections import defaultdict

def analyse_set(fasta_path, set_name):
    records = list(SeqIO.parse(fasta_path, "fasta"))

    total     = len(records)
    lengths   = [len(r.seq) for r in records]
    avg_len   = sum(lengths) / total if total > 0 else 0

    u_counts  = defaultdict(int)
    for r in records:
        seq = str(r.seq).upper()
        u_counts[seq.count('U')] += 1

    print(f"\n{'='*42}")
    print(f" {set_name}")
    print(f"{'='*42}")
    print(f"  Total sequences : {total}")
    print(f"  Avg length      : {avg_len:.0f} aa")
    print(f"  Min length      : {min(lengths)} aa")
    print(f"  Max length      : {max(lengths)} aa")

    if set_name == "POSITIVE SET":
        print(f"\n  Sec (U) count distribution:")
        for n_u, count in sorted(u_counts.items()):
            bar = '█' * count
            print(f"    {n_u} Sec residue(s): {count:>3} proteins  {bar}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--positive', default='data/processed/positive_set.fa')
    parser.add_argument('--negative', default='data/processed/negative_set.fa')
    args = parser.parse_args()
    os.makedirs("data/processed", exist_ok=True)

    analyse_set(args.positive, "POSITIVE SET (Selenoproteins)")
    analyse_set(args.negative, "NEGATIVE SET (Cys-homologs)")

    print(f"\n{'='*42}")
    print(" Reference sets ready for database building")
    print(f"{'='*42}\n")

if __name__ == '__main__':
    main()