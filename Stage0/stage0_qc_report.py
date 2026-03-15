
# stage0_qc_report.py
# Verifies all Stage 0 outputs are present and valid before Stage 1 begins

import argparse
import os
import sys
from datetime import datetime
from Bio import SeqIO

# ─────────────────────────────────────────────
# 1. FILE CHECKS
# ─────────────────────────────────────────────

def check_file(path, label, required=True):
    """
    Check a file exists and is non-empty.
    Returns (status_ok, size_mb)
    """
    exists   = os.path.exists(path)
    size_mb  = os.path.getsize(path) / 1e6 if exists else 0
    nonempty = size_mb > 0

    if exists and nonempty:
        ok     = True
    elif exists and not nonempty:
        
        ok     = False
    else:
        
        ok     = not required   # only fail if required

    print(f"  {label:<45} {size_mb:>8.2f} MB")
    return ok

def check_directory(path, label, min_files=1):
    """
    Check a directory exists and contains at least min_files files.
    """
    exists     = os.path.isdir(path)
    file_count = len(os.listdir(path)) if exists else 0
    ok         = exists and file_count >= min_files

    
    print(f" {label:<45} {file_count:>5} file(s)")
    return ok

# ─────────────────────────────────────────────
# 2. FASTA CONTENT CHECKS
# ─────────────────────────────────────────────

def check_positive_set(fasta_path):
    """
    Verify positive set:
      - Has sequences
      - All contain U (selenocysteine)
      - Show U count distribution
    """
    print("\n--- Positive Set (Selenoproteins) ---")

    if not os.path.exists(fasta_path):
        print("  ❌ File missing — cannot check contents")
        return False

    records   = list(SeqIO.parse(fasta_path, "fasta"))
    total     = len(records)
    with_U    = [r for r in records if 'U' in str(r.seq).upper()]
    without_U = [r for r in records if 'U' not in str(r.seq).upper()]

    print(f"  Total sequences    : {total}")
    print(f"  With U (Sec)       : {len(with_U)}")

    if without_U:
        print(f"  Without U          : {len(without_U)}  (check these)")
        for r in without_U:
            print(f"    - {r.id}")

    # U count distribution
    from collections import Counter
    u_dist = Counter(str(r.seq).upper().count('U') for r in records)
    print(f"  Sec per protein    : ", end="")
    for n, count in sorted(u_dist.items()):
        print(f"{n}×U → {count} proteins  ", end="")
    print()

    return total > 0 and len(without_U) == 0

def check_negative_set(fasta_path):
    """
    Verify negative set:
      - Has sequences
      - None contain U (these are Cys-homologs)
    """
    print("\n--- Negative Set (Cys-homologs) ---")

    if not os.path.exists(fasta_path):
        print("  ❌ File missing — cannot check contents")
        return False

    records  = list(SeqIO.parse(fasta_path, "fasta"))
    total    = len(records)
    with_U   = [r for r in records if 'U' in str(r.seq).upper()]

    print(f"  Total sequences    : {total}")
    print(f"  Without U (clean)  : {total - len(with_U)}")

    if with_U:
        print(f"  WITH U (problem!)  : {len(with_U)}  (these are selenoproteins, not Cys-homologs)")
        for r in with_U:
            print(f"    - {r.id}")
        return False

    return total > 0

def check_ground_truth(gff3_path):
    """
    Verify ground truth GFF3:
      - Has content
      - Contains gene features
      - Not the full annotation (should be small)
    """
    print("\n--- Ground Truth GFF3 ---")

    if not os.path.exists(gff3_path):
        print("  File missing")
        return False

    gene_count = 0
    feat_count = 0

    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            feat_count += 1
            if parts[2] == 'gene':
                gene_count += 1

    print(f"  Gene features      : {gene_count}")
    print(f"  Total features     : {feat_count}")

    if gene_count == 0:
        print("   WARNING: No genes found — check extract_ground_truth.py output")
        return False

    if feat_count > 10000:
        print("  WARNING: Very large GFF3 — may not be filtered correctly")

    return True

def check_genome(fasta_path, index_path):
    """
    Verify genome FASTA and its .fai index.
    """
    print("\n--- Genome ---")

    if not os.path.exists(fasta_path):
        print(" Genome FASTA missing")
        return False

    # Basic genome stats from .fai (fast — no need to parse whole FASTA)
    if os.path.exists(index_path):
        with open(index_path) as f:
            lines = f.readlines()
        print(f"  Sequences indexed  : {len(lines)}")
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom  = parts[0]
                length = int(parts[1])
                print(f"    {chrom:<6} : {length:,} bp  ({length/1e6:.1f} Mb)")
        print(f"  .fai index         : found")
    else:
        print("  .fai index missing — run: samtools faidx chr22.fa")
        return False

    return True

# ─────────────────────────────────────────────
# 3. DATABASE CHECKS
# ─────────────────────────────────────────────

def check_databases(diamond_dir, mmseqs_dir):
    """
    Verify both alignment databases were built correctly.
    """
    print("\n--- Alignment Databases ---")

    # DIAMOND — should have a .dmnd file
    dmnd_files = [f for f in os.listdir(diamond_dir)
                  if f.endswith('.dmnd')] if os.path.isdir(diamond_dir) else []
    if dmnd_files:
        dmnd_path = os.path.join(diamond_dir, dmnd_files[0])
        dmnd_size = os.path.getsize(dmnd_path) / 1024
        print(f"   DIAMOND DB       : {dmnd_files[0]}  ({dmnd_size:.1f} KB)")
        diamond_ok = True
    else:
        print("   DIAMOND DB       : not found")
        diamond_ok = False

    # MMseqs2 — should have selenoDB.index
    mmseqs_index = os.path.join(mmseqs_dir, 'selenoDB.index')
    if os.path.exists(mmseqs_index):
        mmseqs_files = [f for f in os.listdir(mmseqs_dir)
                        if f.startswith('selenoDB')]
        print(f"   MMseqs2 DB       : selenoDB  ({len(mmseqs_files)} files)")
        mmseqs_ok = True
    else:
        print("   MMseqs2 DB       : selenoDB.index not found")
        mmseqs_ok = False

    return diamond_ok and mmseqs_ok

# ─────────────────────────────────────────────
# 4. FINAL REPORT
# ─────────────────────────────────────────────

def print_final_report(checks):
    """
    Print overall pass/fail and what to fix if anything failed.
    """
    all_passed = all(checks.values())

    print("\n==========================================")
    print(" Stage 0 QC Summary")
    print("==========================================")

    for label, passed in checks.items():
        print(f" {label}")

    print("==========================================")

    if all_passed:
        print(" ALL CHECKS PASSED — Stage 1 ready")
    else:
        print(" SOME CHECKS FAILED — fix before running Stage 1")
        print("\n  What to fix:")
        if not checks.get("Genome + index"):
            print("    → Re-run: python bin/prepare_genome.py")
        if not checks.get("Positive set"):
            print("    → Re-run: python bin/download_positive_set.py")
        if not checks.get("Negative set"):
            print("    → Re-run: python bin/download_negative_set.py")
        if not checks.get("Ground truth GFF3"):
            print("    → Re-run: python bin/extract_ground_truth.py")
        if not checks.get("Alignment databases"):
            print("    → Re-run: python bin/prepare_db.py")

    print("==========================================\n")
    return all_passed

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="QC check for all Stage 0 outputs"
    )
    parser.add_argument('--genome',       default='data/raw/chr22.fa')
    parser.add_argument('--genome-index',       default='data/raw/chr22_clean.fa')
    parser.add_argument('--positive-set', default='data/processed/positive_set.fa')
    parser.add_argument('--negative-set', default='data/processed/negative_set.fa')
    parser.add_argument('--ground-truth', default='data/ground_truth/ground_truth_selenoproteins.gff3')
    parser.add_argument('--diamond-dir',  default='databases/diamond_db')
    parser.add_argument('--mmseqs-dir',   default='databases/mmseqs_db')
    args = parser.parse_args()

    print("\n==========================================")
    print(f" Stage 0 QC Report")
    print(f" {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("==========================================")

    print("\n--- File Presence ---")
    check_file(args.genome,                    "Genome (chr22.fa)")
    check_file(args.genome_index,           "Genome index (chr22.fa.fai)")
    check_file(args.positive_set,              "Positive set (positive_set.fa)")
    check_file(args.negative_set,              "Negative set (negative_set.fa)")
    check_file(args.ground_truth,              "Ground truth GFF3")
    check_directory(args.diamond_dir,          "DIAMOND database dir")
    check_directory(args.mmseqs_dir,           "MMseqs2 database dir")

    # Detailed content checks
    genome_ok   = check_genome(args.genome, args.genome_index + '.fai')
    positive_ok = check_positive_set(args.positive_set)
    negative_ok = check_negative_set(args.negative_set)
    truth_ok    = check_ground_truth(args.ground_truth)
    db_ok       = check_databases(args.diamond_dir, args.mmseqs_dir)

    # Final report
    checks = {
        "Genome + index"      : genome_ok,
        "Positive set"        : positive_ok,
        "Negative set"        : negative_ok,
        "Ground truth GFF3"   : truth_ok,
        "Alignment databases" : db_ok,
    }

    passed = print_final_report(checks)
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()