import argparse
import subprocess
import shutil
import sys
import os
import time
from Bio import SeqIO

# ─────────────────────────────────────────────
# 1. CHECK TOOLS ARE INSTALLED
# ─────────────────────────────────────────────

def check_tools():
    """
    Verify DIAMOND and MMseqs2 are available before attempting
    to build databases. Fail early with a clear message.
    """
    print("\n--- Tool Check ---")

    tools = {
        'diamond': 'DIAMOND',
        'mmseqs':  'MMseqs2',
    }

    all_found = True
    for binary, label in tools.items():
        path = shutil.which(binary)
        if path:
            # Get version
            try:
                result = subprocess.run(
                    [binary, 'version'],
                    capture_output=True, text=True
                )
                version = result.stdout.strip().split('\n')[0]
            except Exception:
                version = 'version unknown'
            print(f"  {label:<10} found at {path}  ({version})")
        else:
            print(f"  {label:<10} NOT FOUND install with: conda install -c bioconda {binary}")
            all_found = False

    if not all_found:
        print("\n  ERROR: Missing required tools. Cannot build databases.")
        sys.exit(1)

# ─────────────────────────────────────────────
# 2. VALIDATE INPUT FASTA
# ─────────────────────────────────────────────

def validate_input(fasta_path):
    """
    Check positive_set.fa is valid and all sequences contain U.
    This is a final safety check — validate_and_filter.py should
    have already ensured this, but we verify again before building.
    """
    print("\n--- Input Validation ---")

    if not os.path.exists(fasta_path):
        print(f"  ERROR: File not found: {fasta_path}")
        sys.exit(1)

    records    = list(SeqIO.parse(fasta_path, "fasta"))
    total      = len(records)
    with_U     = sum(1 for r in records if 'U' in str(r.seq).upper())
    without_U  = total - with_U

    print(f"  Total sequences : {total}")
    print(f"  Have U (Sec)    : {with_U}")
    print(f"  Missing U       : {without_U}")

    if total == 0:
        print("  ERROR: positive_set.fa is empty.")
        sys.exit(1)

    if without_U > 0:
        print(f"  WARNING: {without_U} sequences have no U — these may be annotation errors.")
        print("  Proceeding, but check your positive set.")

    return total

# ─────────────────────────────────────────────
# 3. BUILD DIAMOND DATABASE
# ─────────────────────────────────────────────

def build_diamond_db(input_fa, output_dir):
    """
    Runs: diamond makedb --in positive_set.fa -d seleno
    Output: databases/diamond_db/seleno.dmnd
    """
    os.makedirs(output_dir, exist_ok=True)
    db_path = os.path.join(output_dir, "seleno")

    print("\n--- Building DIAMOND Database ---")
    print(f"  Input  : {input_fa}")
    print(f"  Output : {db_path}.dmnd")

    start = time.time()

    result = subprocess.run(
        ['diamond', 'makedb', '--in', input_fa, '-d', db_path],
        capture_output=True,
        text=True
    )

    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"  ERROR: DIAMOND makedb failed:\n{result.stderr}")
        sys.exit(1)

    # Verify output exists and is non-zero
    dmnd_file = db_path + '.dmnd'
    if not os.path.exists(dmnd_file) or os.path.getsize(dmnd_file) == 0:
        print("  ERROR: DIAMOND database file is missing or empty.")
        sys.exit(1)

    size_kb = os.path.getsize(dmnd_file) / 1024
    print(f" Done in {elapsed:.1f}s — size: {size_kb:.1f} KB")
    return dmnd_file

# ─────────────────────────────────────────────
# 4. BUILD MMseqs2 DATABASE
# ─────────────────────────────────────────────

def build_mmseqs_db(input_fa, output_dir):
    """
    Runs: mmseqs createdb positive_set.fa selenoDB
    Output: databases/mmseqs_db/selenoDB (multiple files)
    """
    os.makedirs(output_dir, exist_ok=True)
    db_path = os.path.join(output_dir, "selenoDB")

    print("\n--- Building MMseqs2 Database ---")
    print(f"  Input  : {input_fa}")
    print(f"  Output : {db_path}")

    start = time.time()

    result = subprocess.run(
        ['mmseqs', 'createdb', input_fa, db_path],
        capture_output=True,
        text=True
    )

    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"  ERROR: MMseqs2 createdb failed:\n{result.stderr}")
        sys.exit(1)

    # MMseqs2 creates multiple files — check the main index exists
    index_file = db_path + '.index'
    if not os.path.exists(index_file) or os.path.getsize(index_file) == 0:
        print("  ERROR: MMseqs2 database index is missing or empty.")
        sys.exit(1)

    # Count files created
    db_files = [f for f in os.listdir(output_dir) if f.startswith('selenoDB')]
    print(f" Done in {elapsed:.1f}s — {len(db_files)} files created")
    return db_path

# ─────────────────────────────────────────────
# 5. FINAL SUMMARY
# ─────────────────────────────────────────────

def print_summary(n_sequences, diamond_db, mmseqs_db):
    print("\n==========================================")
    print(" prepare_db.py — Complete")
    print("==========================================")
    print(f"  Sequences indexed : {n_sequences}")
    print(f"  DIAMOND DB        : {diamond_db}")
    print(f"  MMseqs2 DB        : {mmseqs_db}")
    print("\n  Both databases are ready for Stage 1 alignment.")
    print("==========================================\n")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build DIAMOND and MMseqs2 databases from positive selenoprotein set"
    )
    parser.add_argument('--selenoproteins', required=True,
                        help='Filtered positive set FASTA (positive_set.fa)')
    parser.add_argument('--diamond-dir',    default='databases/diamond_db',
                        help='Output directory for DIAMOND database')
    parser.add_argument('--mmseqs-dir',     default='databases/mmseqs_db',
                        help='Output directory for MMseqs2 database')
    args = parser.parse_args()

    print("\n==========================================")
    print(" prepare_db.py")
    print(" Building Alignment Databases")
    print("==========================================")

    # Step 1 — Check tools
    check_tools()

    # Step 2 — Validate input
    n_sequences = validate_input(args.selenoproteins)

    # Step 3 — Build DIAMOND DB
    diamond_db = build_diamond_db(args.selenoproteins, args.diamond_dir)

    # Step 4 — Build MMseqs2 DB
    mmseqs_db = build_mmseqs_db(args.selenoproteins, args.mmseqs_dir)

    # Step 5 — Summary
    print_summary(n_sequences, diamond_db, mmseqs_db)


if __name__ == '__main__':
    main()