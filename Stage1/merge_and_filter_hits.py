
# merge_and_filter_hits.py
# Stage 1: Merge DIAMOND + MMseqs2 hits, apply quality filters, deduplicate

import argparse
import os
import sys
import pandas as pd

DEFAULT_MIN_IDENTITY = 50.0
DEFAULT_MIN_COVERAGE = 5.0
DEFAULT_MAX_EVALUE   = 1e-5
DEFAULT_MIN_BITSCORE = 60.0   # eliminates the MMseqs2 noise (all scored <54)
OVERLAP_TOLERANCE    = 50     # bp

# ─────────────────────────────────────────────
# 1. LOAD HITS
# ─────────────────────────────────────────────

def load_hits(diamond_tsv, mmseqs_tsv):
    print("\n--- Loading Hits ---")

    dfs = []

    for path, label in [(diamond_tsv, 'DIAMOND'), (mmseqs_tsv, 'MMseqs2')]:
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found — skipping {label}")
            continue

        df = pd.read_csv(path, sep='\t')

        if df.empty:
            print(f"  {label:<10} : 0 hits")
            continue

        dfs.append(df)
        print(f"  {label:<10} : {len(df)} hits loaded")

    if not dfs:
        print("  ERROR: No hits from either tool.")
        sys.exit(1)

    combined = pd.concat(dfs, ignore_index=True)
    print(f"  Combined   : {len(combined)} total hits")
    return combined

# ─────────────────────────────────────────────
# 2. NORMALISE PROTEIN IDs
# ─────────────────────────────────────────────

def normalise_protein_ids(df):
    """
    DIAMOND uses full UniProt format: sp|Q16881|TRXR1_HUMAN
    MMseqs2 uses accession only:      Q16881

    Normalise both to accession only so cross-tool
    deduplication works correctly.

    sp|Q16881|TRXR1_HUMAN  →  Q16881
    Q16881                 →  Q16881  (unchanged)
    """
    print(f"\n--- Normalising Protein IDs ---")

    def extract_accession(pid):
        pid = str(pid)
        if pid.startswith('sp|') or pid.startswith('tr|'):
            parts = pid.split('|')
            if len(parts) >= 2:
                return parts[1]
        return pid

    df = df.copy()
    df['target_protein_raw'] = df['target_protein']          # keep original
    df['target_protein']     = df['target_protein'].apply(extract_accession)

    # before/after for first few
    examples = df[['target_protein_raw', 'target_protein']].drop_duplicates().head(3)
    for _, row in examples.iterrows():
        if row['target_protein_raw'] != row['target_protein']:
            print(f"  {row['target_protein_raw'][:40]:<40} → {row['target_protein']}")

    print(f"  Unique proteins after normalisation: {df['target_protein'].nunique()}")
    return df

# ─────────────────────────────────────────────
# 3. APPLY QUALITY FILTERS
# ─────────────────────────────────────────────

def apply_filters(df, min_identity, min_coverage,
                  max_evalue, min_bitscore):
    """
    Four quality filters applied in order:
      1. Identity  >= min_identity   (default 50%)
      2. Coverage  >= min_coverage   (default 5% — exon fragments are normal)
      3. E-value   <= max_evalue     (default 1e-5)
      4. Bitscore  >= min_bitscore   (default 60 — eliminates repeat noise)

    Note on coverage:
      Genomic DNA → protein blastx finds individual exons.
      Each exon covers only a small fraction of the full protein.
      Real hits: SELM exon = 36%, SELENOO exon = 28%, TXNRD2 exon = 13%
      Setting coverage >= 50% would discard all real selenoprotein hits.
      Use bitscore as the primary quality signal instead.
    """
    print(f"\n--- Applying Quality Filters ---")
    print(f"  Min identity : {min_identity}%")
    print(f"  Min coverage : {min_coverage}%")
    print(f"  Max e-value  : {max_evalue:.0e}")
    print(f"  Min bitscore : {min_bitscore}")

    before = len(df)

    after_identity = df[df['identity'] >= min_identity]
    after_coverage = after_identity[after_identity['coverage'] >= min_coverage]
    after_evalue   = after_coverage[after_coverage['evalue'] <= max_evalue]
    after_bitscore = after_evalue[after_evalue['bitscore'] >= min_bitscore]

    print(f"\n  Before filters        : {before}")
    print(f"  After identity filter : {len(after_identity)}"
          f"  (-{before - len(after_identity)})")
    print(f"  After coverage filter : {len(after_coverage)}"
          f"  (-{len(after_identity) - len(after_coverage)})")
    print(f"  After e-value filter  : {len(after_evalue)}"
          f"  (-{len(after_coverage) - len(after_evalue)})")
    print(f"  After bitscore filter : {len(after_bitscore)}"
          f"  (-{len(after_evalue) - len(after_bitscore)})  ← removes repeat noise")

    return after_bitscore.copy()

# ─────────────────────────────────────────────
# 4. DEDUPLICATE OVERLAPPING HITS  (fixed O(n log n))
# ─────────────────────────────────────────────

def deduplicate_hits(df):
    """
    Mark lower-scoring overlapping hits as redundant.

    Fixed issues vs original:
      1. Uses normalised protein IDs for cross-tool matching
      2. Groups by (chrom, protein) before overlap check → O(n log n)
         instead of O(n²) — critical for large hit sets
    """
    print(f"\n--- Deduplicating Overlapping Hits ---")

    df         = df.copy()
    df['overlap_group'] = -1
    df['redundant']     = False

    group_id  = 0
    processed = set()

    # Group by chrom + protein first — only compare within same group
    # This reduces comparisons from O(n²) to O(k²) per group
    grouped = df.groupby(['query_id', 'target_protein'])

    for (chrom, protein), group in grouped:
        indices = list(group.index)

        if len(indices) == 1:
            continue   # nothing to deduplicate

        # Sort by bitscore descending — best hit first
        group_sorted = group.sort_values('bitscore', ascending=False)
        sorted_idx   = list(group_sorted.index)

        local_processed = set()

        for i, idx_i in enumerate(sorted_idx):
            if idx_i in processed or idx_i in local_processed:
                continue

            row_i       = df.loc[idx_i]
            overlapping = [idx_i]

            for idx_j in sorted_idx[i+1:]:
                if idx_j in processed or idx_j in local_processed:
                    continue

                row_j = df.loc[idx_j]

                overlap = (
                    row_i['query_start'] <= row_j['query_end'] + OVERLAP_TOLERANCE
                    and
                    row_j['query_start'] <= row_i['query_end'] + OVERLAP_TOLERANCE
                )

                if overlap:
                    overlapping.append(idx_j)

            if len(overlapping) > 1:
                df.loc[overlapping, 'overlap_group'] = group_id
                # First in sorted list = highest bitscore = keep
                for idx in overlapping[1:]:
                    df.loc[idx, 'redundant'] = True
                group_id += 1

            local_processed.update(overlapping)

        processed.update(sorted_idx)

    n_redundant = df['redundant'].sum()
    n_unique    = len(df) - n_redundant

    print(f"  Overlap groups found  : {group_id}")
    print(f"  Redundant hits flagged: {n_redundant}")
    print(f"  Unique hits remaining : {n_unique}")

    return df

# ─────────────────────────────────────────────
# 5. TOOL COMPARISON SUMMARY
# ─────────────────────────────────────────────

def compare_tools(df):
    print(f"\n--- Tool Comparison ---")

    for tool in ['DIAMOND', 'MMseqs2']:
        subset = df[df['source'] == tool]
        if subset.empty:
            print(f"  {tool:<10} : no hits after filtering")
            continue

        print(f"\n  {tool}:")
        print(f"    Hits              : {len(subset)}")
        print(f"    Unique proteins   : {subset['target_protein'].nunique()}")
        print(f"    Mean identity     : {subset['identity'].mean():.1f}%")
        print(f"    Mean coverage     : {subset['coverage'].mean():.1f}%")
        print(f"    Mean bitscore     : {subset['bitscore'].mean():.1f}")

    # Use normalised protein IDs for overlap comparison
    diamond_proteins = set(df[df['source'] == 'DIAMOND']['target_protein'])
    mmseqs_proteins  = set(df[df['source'] == 'MMseqs2']['target_protein'])

    both     = diamond_proteins & mmseqs_proteins
    only_dia = diamond_proteins - mmseqs_proteins
    only_mm  = mmseqs_proteins  - diamond_proteins

    print(f"\n  Protein overlap (using normalised IDs):")
    print(f"    Found by both tools  : {len(both)}")
    print(f"    DIAMOND only         : {len(only_dia)}")
    print(f"    MMseqs2 only         : {len(only_mm)}")

    if only_dia:
        print(f"\n    Proteins unique to DIAMOND:")
        for p in sorted(only_dia):
            print(f"      {p}")
    if only_mm:
        print(f"\n    Proteins unique to MMseqs2:")
        for p in sorted(only_mm):
            print(f"      {p}")

# ─────────────────────────────────────────────
# 6. WRITE OUTPUTS
# ─────────────────────────────────────────────

def write_outputs(df, output_tsv, output_summary):
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n  Merged hits saved    : {output_tsv}")

    summary_rows = []
    for tool in ['DIAMOND', 'MMseqs2', 'combined']:
        if tool == 'combined':
            subset = df[df['redundant'] == False]
            label  = 'combined (unique)'
        else:
            subset = df[df['source'] == tool]
            label  = tool

        if subset.empty:
            continue

        summary_rows.append({
            'tool'            : label,
            'total_hits'      : len(subset),
            'unique_proteins' : subset['target_protein'].nunique(),
            'mean_identity'   : round(subset['identity'].mean(), 2),
            'mean_coverage'   : round(subset['coverage'].mean(), 2),
            'mean_bitscore'   : round(subset['bitscore'].mean(), 2),
        })

    pd.DataFrame(summary_rows).to_csv(output_summary, sep='\t', index=False)
    print(f"  Summary saved        : {output_summary}")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 1: Merge DIAMOND + MMseqs2 hits and apply filters"
    )
    parser.add_argument('--diamond-hits',    default='results/stage1/hits_diamond.tsv')
    parser.add_argument('--mmseqs-hits',     default='results/stage1/hits_mmseqs2.tsv')
    parser.add_argument('--output-tsv',      default='results/stage1/merged_hits.tsv')
    parser.add_argument('--output-summary',  default='results/stage1/merge_summary.tsv')
    parser.add_argument('--min-identity',    type=float, default=DEFAULT_MIN_IDENTITY)
    parser.add_argument('--min-coverage',    type=float, default=DEFAULT_MIN_COVERAGE)
    parser.add_argument('--max-evalue',      type=float, default=DEFAULT_MAX_EVALUE)
    parser.add_argument('--min-bitscore',    type=float, default=DEFAULT_MIN_BITSCORE)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    print("\n==========================================")
    print(" Stage 1: Merge and Filter Hits")
    print("==========================================")

    df = load_hits(args.diamond_hits, args.mmseqs_hits)
    df = normalise_protein_ids(df)
    df = apply_filters(df, args.min_identity, args.min_coverage,
                       args.max_evalue, args.min_bitscore)

    if df.empty:
        print("\n  WARNING: No hits passed filters.")
        sys.exit(0)

    df = deduplicate_hits(df)
    compare_tools(df)
    write_outputs(df, args.output_tsv, args.output_summary)

    unique_hits = df[df['redundant'] == False]

    print(f"\n==========================================")
    print(f" Stage 1 Merge — Complete")
    print(f"==========================================")
    print(f"  Total hits (all)     : {len(df)}")
    print(f"  Unique hits          : {len(unique_hits)}")
    print(f"  Proteins covered     : {df['target_protein'].nunique()}")
    print(f"  Ready for Stage 2    : {args.output_tsv}")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()