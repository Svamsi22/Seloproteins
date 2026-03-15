import argparse
import subprocess
import shutil
import sys
import os
from Bio import SeqIO

# ─────────────────────────────────────────────
# 1. CHECK MASKING
# ─────────────────────────────────────────────

def check_masking(fasta_path):
    """
    Reads first 10,000 bases of genome to detect masking type.

    Soft-masked  → repeats in lowercase  (best for alignment)
    Hard-masked  → repeats replaced with N
    Unmasked     → all uppercase, no N's
    """
    sample_seq = ""
    for record in SeqIO.parse(fasta_path, "fasta"):
        sample_seq = str(record.seq[:10000])
        break

    lowercase_count = sum(1 for c in sample_seq if c.islower())
    n_count         = sample_seq.upper().count('N')

    print("\n--- Masking Check ---")

    if lowercase_count > 0:
        mask_type = "SOFT-MASKED"
        status    = "OK — recommended for alignment"
    elif n_count > 100:
        mask_type = "HARD-MASKED"
        status    = "WARNING — may miss selenoproteins in repeat regions"
    else:
        mask_type = "UNMASKED"
        status    = "WARNING — may produce false positives in repeat regions"

    print(f"  Type   : {mask_type}")
    print(f"  Status : {status}")

    return mask_type

# ─────────────────────────────────────────────
# 2. CHECK CHROMOSOME NAMING
# ─────────────────────────────────────────────

def check_genome_naming(fasta_path):
    """
    Detects whether genome uses UCSC (chr22) or Ensembl (22) naming.
    Returns the first chromosome name found as an example.
    """
    for record in SeqIO.parse(fasta_path, "fasta"):
        first_id       = record.id
        has_chr_prefix = first_id.startswith('chr')

        print("\n--- Chromosome Naming Check (Genome) ---")
        print(f"  Example chrom : {first_id}")
        print(f"  Style         : "
              f"{'UCSC (chr prefix)' if has_chr_prefix else 'Ensembl (no prefix)'}")

        return has_chr_prefix

def check_gff3_naming(gff3_path):
    """
    Detects whether GFF3 uses UCSC or Ensembl chromosome naming.
    """
    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                chrom          = parts[0]
                has_chr_prefix = chrom.startswith('chr')

                print("\n--- Chromosome Naming Check (GFF3) ---")
                print(f"  Example chrom : {chrom}")
                print(f"  Style         : "
                      f"{'UCSC (chr prefix)' if has_chr_prefix else 'Ensembl (no prefix)'}")

                return has_chr_prefix

    return False

# ─────────────────────────────────────────────
# 3. STANDARDISE NAMING IF MISMATCH
# ─────────────────────────────────────────────

def standardise_chromosome_names(fasta_path, output_path,
                                  add_prefix=False,
                                  remove_prefix=False):
    """
    Renames chromosomes in FASTA to match GFF3 convention.
    Only runs if mismatch is detected between genome and GFF3.

    add_prefix    → Ensembl genome + UCSC GFF3  → add 'chr'
    remove_prefix → UCSC genome + Ensembl GFF3  → remove 'chr'
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))

    for record in records:
        if add_prefix and not record.id.startswith('chr'):
            record.id          = 'chr' + record.id
            record.description = record.id

        elif remove_prefix and record.id.startswith('chr'):
            record.id          = record.id[3:]
            record.description = record.id

    SeqIO.write(records, output_path, "fasta")
    print(f"\n  Chromosome names updated → {output_path}")

# ─────────────────────────────────────────────
# 4. INDEX GENOME
# ─────────────────────────────────────────────

def index_genome(fasta_path):
    """
    Runs samtools faidx to create .fai index.
    Required by BEDTools in Stage 2 for region extraction.
    """
    print("\n--- Genome Indexing ---")

    # Check samtools is available
    if not shutil.which('samtools'):
        print("  ERROR: samtools not found. Please install samtools.")
        sys.exit(1)

    result = subprocess.run(
        ['samtools', 'faidx', fasta_path],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"  ERROR: samtools faidx failed:\n{result.stderr}")
        sys.exit(1)

    index_path = fasta_path + '.fai'
    index_size = os.path.getsize(index_path)

    print(f"  Index created : {index_path}")
    print(f"  Index size    : {index_size:,} bytes")

# ─────────────────────────────────────────────
# 5. GENOME STATS
# ─────────────────────────────────────────────

def print_genome_stats(fasta_path):
    """
    Prints basic stats about the genome file.
    """
    total_bases = 0
    n_count     = 0
    sequences   = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq          = str(record.seq).upper()
        total_bases += len(seq)
        n_count     += seq.count('N')
        sequences   += 1

    n_percent = (n_count / total_bases * 100) if total_bases > 0 else 0

    print("\n--- Genome Stats ---")
    print(f"  Sequences     : {sequences}")
    print(f"  Total bases   : {total_bases:,} bp")
    print(f"  Size          : {total_bases / 1e6:.1f} Mb")
    print(f"  N bases       : {n_count:,} ({n_percent:.1f}%)")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Validate, standardise, and index genome for pipeline"
    )
    parser.add_argument('--genome',  default='data/raw/chr22.fa',
                        help='Input genome FASTA (chr22.fa)')
    parser.add_argument('--gff3',    default='data/ground_truth/ensembl_chr22_full.gff3',
                        help='Ground truth GFF3 to check naming against')
    parser.add_argument('--output',  default='data/raw/chr22_clean.fa',
                        help='Output genome FASTA (standardised)')
    args = parser.parse_args()

    print("\n==========================================")
    print(" prepare_genome.py")
    print("==========================================")

    # Step 1 — Check masking
    check_masking(args.genome)

    # Step 2 — Check naming convention in genome vs GFF3
    genome_has_chr = check_genome_naming(args.genome)
    gff3_has_chr   = check_gff3_naming(args.gff3)

    # Step 3 — Fix naming if mismatch
    print("\n--- Naming Standardisation ---")
    if genome_has_chr == gff3_has_chr:
        print("  Naming is consistent — no changes needed")
        shutil.copy(args.genome, args.output)
        print(f"  Copied as-is → {args.output}")
    else:
        print("  Mismatch detected — standardising...")
        standardise_chromosome_names(
            args.genome,
            args.output,
            add_prefix    = (gff3_has_chr and not genome_has_chr),
            remove_prefix = (genome_has_chr and not gff3_has_chr)
        )

    # Step 4 — Index genome
    index_genome(args.output)

    # Step 5 — Print stats
    print_genome_stats(args.output)

    print("\n==========================================")
    print(" prepare_genome.py — Complete")
    print("==========================================")
    print(f"  Output genome : {args.output}")
    print(f"  Genome index  : {args.output}.fai")
    print("==========================================\n")

if __name__ == '__main__':
    main()