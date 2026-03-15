
# extract_sequences.py
# Stage 2: Extract FASTA sequences for each extended region
# Uses bedtools getfasta to pull sequences from the genome

import argparse
import os
import sys
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO

# ─────────────────────────────────────────────
# WHY THIS SCRIPT EXISTS
# ─────────────────────────────────────────────
# Stage 3 needs actual DNA sequences to:
#   1. Check for TGA (UGA) codons at Sec positions
#   2. Run SECIS element detection on downstream region
#
# This script pulls those sequences from the genome
# using the extended coordinates from extend_regions.py
#
# Output is a FASTA where each record = one candidate region
# The sequence includes:
#   [upstream_buffer][gene][downstream_buffer_for_SECIS]

# ─────────────────────────────────────────────
# 1. CHECK BEDTOOLS
# ─────────────────────────────────────────────

def check_bedtools():
    if not shutil.which('bedtools'):
        print("  ERROR: BEDTools not found.")
        print("  Install with: conda install -c bioconda bedtools")
        sys.exit(1)

    result  = subprocess.run(
        ['bedtools', '--version'],
        capture_output=True, text=True
    )
    print(f"  BEDTools: {result.stdout.strip()}")

# ─────────────────────────────────────────────
# 2. RUN BEDTOOLS GETFASTA
# ─────────────────────────────────────────────

def run_getfasta(genome_fa, extended_bed, output_fa):
    """
    bedtools getfasta extracts sequences from genome FASTA
    using coordinates in the BED file.

    Flags used:
      -fi   : input genome FASTA
      -bed  : BED file with coordinates
      -fo   : output FASTA
      -s    : strand-aware extraction
              reverse-complements - strand sequences
              so all output sequences read 5'→3'
      -name : use BED name field as FASTA header
              (our name field = ext_region_id)
    """
    print(f"\n--- Running bedtools getfasta ---")
    print(f"  Genome  : {genome_fa}")
    print(f"  BED     : {extended_bed}")
    print(f"  Output  : {output_fa}")

    cmd = [
        'bedtools', 'getfasta',
        '-fi',   genome_fa,
        '-bed',  extended_bed,
        '-fo',   output_fa,
        '-s',                   # strand-aware
        '-name',                # use name col as header
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: bedtools getfasta failed:\n{result.stderr}")
        sys.exit(1)

    # bedtools sometimes prints warnings to stderr even on success
    if result.stderr:
        print(f"  WARNING from bedtools: {result.stderr.strip()}")

    print(f"  Done.")

# ─────────────────────────────────────────────
# 3. VALIDATE FASTA OUTPUT
# ─────────────────────────────────────────────

def validate_fasta(output_fa):
    """
    Check the extracted FASTA is valid and non-empty.

    Common failure modes:
      - Empty file: BED coordinates outside chromosome range
      - Wrong chrom names: silent mismatch between BED and genome
      - N-only sequences: region fell entirely in hard-masked area
    """
    print(f"\n--- Validating Extracted FASTA ---")

    if not os.path.exists(output_fa) or os.path.getsize(output_fa) == 0:
        print(f"  ERROR: Output FASTA is empty.")
        print(f"  Check chromosome naming and coordinates.")
        sys.exit(1)

    records    = list(SeqIO.parse(output_fa, "fasta"))
    total      = len(records)
    empty_seqs = [r for r in records if len(r.seq) == 0]
    n_only     = [r for r in records
                  if len(r.seq) > 0 and
                  str(r.seq).upper().replace('N', '') == '']

    print(f"  Sequences extracted  : {total}")

    if empty_seqs:
        print(f"  Empty sequences    : {len(empty_seqs)}")
        for r in empty_seqs:
            print(f"    {r.id}")

    if n_only:
        print(f"  All-N sequences    : {len(n_only)}")
        print(f"     These regions are hard-masked — no signal possible")
        for r in n_only:
            print(f"    {r.id}")

    valid = total - len(empty_seqs) - len(n_only)
    print(f"  Valid sequences      : {valid}")

    return records

# ─────────────────────────────────────────────
# 4. ANNOTATE FASTA WITH METADATA
# ─────────────────────────────────────────────

def annotate_fasta(records, extended_tsv, output_annotated_fa):

    ext_df = pd.read_csv(extended_tsv, sep='\t')

    meta = {}
    for _, row in ext_df.iterrows():
        ext_id = row['ext_region_id']
        meta[ext_id]                = row
        meta[f"{ext_id}::{ext_id}"] = row

    annotated = []
    matched   = 0
    unmatched = 0

    for record in records:
        region_id        = record.id
        region_id_simple = region_id.split('::')[0]

        # ── explicit None check (not 'or') ──────────────────
        row_meta = meta.get(region_id)
        if row_meta is None:
            row_meta = meta.get(region_id_simple)
        # ─────────────────────────────────────────────────────

        if row_meta is not None:
            record.description = (
                f"target={row_meta['target_protein']} "
                f"strand={row_meta['strand']} "
                f"gene={row_meta['chrom']}:{row_meta['gene_start']}-{row_meta['gene_end']} "
                f"identity={row_meta['identity']} "
                f"coverage={row_meta['coverage']} "
                f"source={row_meta['sources']}"
            )
            matched += 1
        else:
            unmatched += 1

        annotated.append(record)

    SeqIO.write(annotated, output_annotated_fa, "fasta")
    print(f"  Matched to metadata  : {matched}")
    print(f"  Unmatched (no meta)  : {unmatched}")
    print(f"  Annotated FASTA saved: {output_annotated_fa}")
    return annotated

# ─────────────────────────────────────────────
# 5. WRITE SEQUENCE STATS TSV
# ─────────────────────────────────────────────

def write_sequence_stats(records, extended_tsv, output_stats):
    print(f"\n--- Writing Sequence Stats ---")

    ext_df = pd.read_csv(extended_tsv, sep='\t')

    meta = {}
    for _, row in ext_df.iterrows():
        ext_id = row['ext_region_id']
        meta[ext_id]                = row
        meta[f"{ext_id}::{ext_id}"] = row

    rows = []
    for record in records:
        seq    = str(record.seq).upper()
        length = len(seq)

        if length == 0:
            continue

        gc    = (seq.count('G') + seq.count('C')) / length * 100
        n_pct = seq.count('N') / length * 100

        region_id        = record.id
        region_id_simple = region_id.split('::')[0]

        # ── explicit None check ──────────────────────────────
        row_meta = meta.get(region_id)
        if row_meta is None:
            row_meta = meta.get(region_id_simple)
        if row_meta is None:
            row_meta = {}
        # ─────────────────────────────────────────────────────

        rows.append({
            'region_id'      : record.id,
            'seq_length'     : length,
            'gc_content'     : round(gc, 2),
            'n_percent'      : round(n_pct, 2),
            'target_protein' : row_meta.get('target_protein', ''),
            'strand'         : row_meta.get('strand', ''),
            'gene_start'     : row_meta.get('gene_start', ''),
            'gene_end'       : row_meta.get('gene_end', ''),
            'sources'        : row_meta.get('sources', ''),
            'identity'       : row_meta.get('identity', ''),
            'coverage'       : row_meta.get('coverage', ''),
        })

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(output_stats, sep='\t', index=False)
    print(f"  Stats saved          : {output_stats}")
    return stats_df

# ─────────────────────────────────────────────
# 6. SUMMARY
# ─────────────────────────────────────────────

def print_summary(stats_df):
    print(f"\n--- Sequence Extraction Summary ---")

    if stats_df.empty:
        print("  No sequences extracted.")
        return

    print(f"  Total sequences      : {len(stats_df)}")
    print(f"  Seq length (min)     : {stats_df['seq_length'].min():,} bp")
    print(f"  Seq length (max)     : {stats_df['seq_length'].max():,} bp")
    print(f"  Seq length (mean)    : {stats_df['seq_length'].mean():,.0f} bp")
    print(f"  Mean GC content      : {stats_df['gc_content'].mean():.1f}%")
    print(f"  Mean N content       : {stats_df['n_percent'].mean():.1f}%")

    print(f"\n  Sequences per source:")
    if 'sources' in stats_df.columns:
        for source, count in stats_df['sources'].value_counts().items():
            print(f"    {source:<20} : {count}")

    print(f"\n  Per-region breakdown:")
    print(f"  {'Region':<40} {'Length':>8}  {'GC%':>6}  {'N%':>5}  {'Target'}")
    print(f"  {'-'*40} {'-'*8}  {'-'*6}  {'-'*5}  {'-'*20}")
    for _, row in stats_df.iterrows():
        print(
            f"  {str(row['region_id']):<40} "
            f"{int(row['seq_length']):>8,}  "
            f"{row['gc_content']:>6.1f}  "
            f"{row['n_percent']:>5.1f}  "
            f"{str(row['target_protein'])}"
        )

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: Extract candidate region sequences from genome"
    )
    parser.add_argument('--genome',
                        default='data/raw/chr22_clean.fa',
                        help='Genome FASTA')
    parser.add_argument('--extended-bed',
                        default='results/stage2/extended_regions.bed',
                        help='Extended BED from extend_regions.py')
    parser.add_argument('--extended-tsv',
                        default='results/stage2/extended_regions.tsv',
                        help='Extended TSV with metadata')
    parser.add_argument('--output-fa',
                        default='results/stage2/candidate_regions_raw.fa',
                        help='Raw extracted FASTA (bedtools output)')
    parser.add_argument('--output-annotated-fa',
                        default='results/stage2/candidate_regions.fa',
                        help='Annotated FASTA with metadata in headers')
    parser.add_argument('--output-stats',
                        default='results/stage2/candidate_regions_stats.tsv',
                        help='Per-sequence stats TSV')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_fa), exist_ok=True)

    print("\n==========================================")
    print(" Stage 2: Extract Sequences")
    print("==========================================")

    # Step 1 — Check BEDTools
    check_bedtools()

    # Step 2 — Run bedtools getfasta
    run_getfasta(args.genome, args.extended_bed, args.output_fa)

    # Step 3 — Validate output
    records = validate_fasta(args.output_fa)

    # Step 4 — Annotate headers with metadata
    annotated = annotate_fasta(records, args.extended_tsv,
                               args.output_annotated_fa)

    # Step 5 — Write sequence stats
    stats_df = write_sequence_stats(annotated, args.extended_tsv,
                                    args.output_stats)

    # Step 6 — Summary
    print_summary(stats_df)

    print(f"\n==========================================")
    print(f" Stage 2 extract_sequences — Complete")
    print(f"==========================================")
    print(f"  Candidate FASTA : {args.output_annotated_fa}")
    print(f"  Sequence stats  : {args.output_stats}")
    print(f"  Ready for Stage 3 validation")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()