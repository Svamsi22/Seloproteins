
# merge_nearby_hits.py
# Stage 2: Merge nearby/overlapping BED hits before region extraction
# Prevents extracting the same gene region multiple times

import argparse
import os
import sys
import shutil
import subprocess
import pandas as pd

# ─────────────────────────────────────────────
# MERGE DISTANCE
# ─────────────────────────────────────────────
# Hits within this distance (bp) are merged into one region.
#
# Why 500bp default?
#   Exons of the same gene are typically split across
#   multiple alignment hits. Gaps between exons in
#   human genes are usually < 500bp for short introns.
#   Larger gaps = likely different genes, don't merge.
#
# Example:
#   Hit 1:  19875520 - 19875680   (exon 1)
#   Hit 2:  19875820 - 19875950   (exon 2, 140bp away)
#   Hit 3:  19876110 - 19876250   (exon 3, 160bp away)
#   Merged: 19875520 - 19876250   (whole gene region)

DEFAULT_MERGE_DISTANCE = 500   # bp

# ─────────────────────────────────────────────
# 1. CHECK BEDTOOLS
# ─────────────────────────────────────────────

def check_bedtools():
    if not shutil.which('bedtools'):
        print("  ERROR: BEDTools not found.")
        sys.exit(1)

    result  = subprocess.run(
        ['bedtools', '--version'],
        capture_output=True, text=True
    )
    version = result.stdout.strip()
    print(f"  BEDTools version: {version}")

# ─────────────────────────────────────────────
# 2. SORT BED FILE
# ─────────────────────────────────────────────

def sort_bed(input_bed, sorted_bed):
    """
    BEDTools merge requires input to be sorted by
    chromosome then start position.
    bedtools sort handles this correctly.
    """
    print(f"\n--- Sorting BED File ---")

    cmd = [
        'bedtools', 'sort',
        '-i', input_bed
    ]

    with open(sorted_bed, 'w') as out:
        result = subprocess.run(cmd, stdout=out, capture_output=False)

    if result.returncode != 0:
        print(f"  ERROR: bedtools sort failed")
        sys.exit(1)

    n_lines = sum(1 for _ in open(sorted_bed))
    print(f"  Sorted BED saved : {sorted_bed}")
    print(f"  Entries          : {n_lines}")

# ─────────────────────────────────────────────
# 3. MERGE NEARBY HITS
# ─────────────────────────────────────────────

def merge_hits(sorted_bed, merged_bed, merge_distance):
    """
    Run bedtools merge to collapse nearby hits.

    -d merge_distance  : merge hits within this many bp
    -s                 : strand-aware merging
                         hits on + and - strand NOT merged together
                         (important — same locus, opposite strands
                          = different genes, must stay separate)
    -c 4,5,6           : columns to carry over from input
                         col 4 = name (metadata)
                         col 5 = score (bitscore)
                         col 6 = strand
    -o collapse,max,distinct : how to handle multiple values per column
                         name   → collapse  (join all names with comma)
                         score  → max       (keep highest bitscore)
                         strand → distinct  (should all be same after -s)
    """
    print(f"\n--- Merging Nearby Hits ---")
    print(f"  Merge distance : {merge_distance} bp")
    print(f"  Strand-aware   : yes")

    cmd = [
        'bedtools', 'merge',
        '-i',  sorted_bed,
        '-d',  str(merge_distance),
        '-s',                          # strand-aware
        '-c',  '4,5,6',               # carry over name, score, strand
        '-o',  'collapse,max,distinct' # how to combine them
    ]

    with open(merged_bed, 'w') as out:
        result = subprocess.run(cmd, stdout=out, capture_output=False)

    if result.returncode != 0:
        print(f"  ERROR: bedtools merge failed")
        sys.exit(1)

    n_lines = sum(1 for _ in open(merged_bed))
    print(f"  Merged BED saved : {merged_bed}")
    print(f"  Regions after merge : {n_lines}")

# ─────────────────────────────────────────────
# 4. PARSE MERGED BED + REPORT
# ─────────────────────────────────────────────

def parse_merged_bed(merged_bed, output_tsv):
    """
    Parse the merged BED into a clean TSV with proper column names.

    After bedtools merge -c 4,5,6 -o collapse,max,distinct:
      col 0: chrom
      col 1: start
      col 2: end
      col 3: names     (comma-separated, all original hit names)
      col 4: max_score (highest bitscore among merged hits)
      col 5: strand

    We also:
      - Parse the encoded name fields back into metadata
      - Calculate region size
      - Count how many hits were merged per region
    """
    print(f"\n--- Parsing Merged Regions ---")

    rows = []
    with open(merged_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            chrom     = parts[0]
            start     = int(parts[1])
            end       = int(parts[2])
            names     = parts[3]          # collapsed names
            max_score = float(parts[4])
            strand    = parts[5]

            # Parse the first name entry to recover metadata
            # Format: query_id|target_protein|source|identity|coverage
            first_name  = names.split(',')[0]
            name_parts  = first_name.split('|')

            if len(name_parts) == 5:
                target_protein = name_parts[1]
                sources        = ','.join(set(
                    n.split('|')[2] for n in names.split(',')
                    if len(n.split('|')) == 5
                ))
                identity       = float(name_parts[3])
                coverage       = float(name_parts[4])
            else:
                target_protein = 'unknown'
                sources        = 'unknown'
                identity       = 0.0
                coverage       = 0.0

            n_hits_merged = len(names.split(','))
            region_size   = end - start

            # Build unique region ID for tracking through pipeline
            region_id = f"{chrom}:{start}-{end}({strand})"

            rows.append({
                'region_id'      : region_id,
                'chrom'          : chrom,
                'start'          : start,
                'end'            : end,
                'strand'         : strand,
                'target_protein' : target_protein,
                'sources'        : sources,
                'identity'       : identity,
                'coverage'       : coverage,
                'max_bitscore'   : max_score,
                'n_hits_merged'  : n_hits_merged,
                'region_size_bp' : region_size,
                'raw_names'      : names,
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_tsv, sep='\t', index=False)

    print(f"  Regions parsed       : {len(df)}")
    print(f"  Saved to             : {output_tsv}")

    return df

# ─────────────────────────────────────────────
# 5. SUMMARY
# ─────────────────────────────────────────────

def print_summary(input_bed, df, merge_distance):
    """
    Show before/after merge counts and what was collapsed.
    """
    n_before = sum(1 for _ in open(input_bed))
    n_after  = len(df)
    n_merged = n_before - n_after

    print(f"\n--- Merge Summary ---")
    print(f"  Hits before merge    : {n_before}")
    print(f"  Regions after merge  : {n_after}")
    print(f"  Hits collapsed       : {n_merged}")
    print(f"  Merge distance used  : {merge_distance} bp")

    if df.empty:
        return

    # Show regions where multiple hits were merged
    multi = df[df['n_hits_merged'] > 1]
    if not multi.empty:
        print(f"\n  Regions with multiple merged hits ({len(multi)}):")
        for _, row in multi.iterrows():
            print(f"    {row['region_id']}")
            print(f"      Protein  : {row['target_protein']}")
            print(f"      Merged   : {row['n_hits_merged']} hits")
            print(f"      Size     : {row['region_size_bp']:,} bp")
            print(f"      Sources  : {row['sources']}")

    print(f"\n  Region size stats:")
    print(f"    Min  : {df['region_size_bp'].min():,} bp")
    print(f"    Max  : {df['region_size_bp'].max():,} bp")
    print(f"    Mean : {df['region_size_bp'].mean():,.0f} bp")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: Merge nearby BED hits before region extraction"
    )
    parser.add_argument('--input-bed',
                        default='results/stage2/hits.bed',
                        help='BED file from hits_to_bed.py')
    parser.add_argument('--sorted-bed',
                        default='results/stage2/hits_sorted.bed',
                        help='Intermediate sorted BED')
    parser.add_argument('--merged-bed',
                        default='results/stage2/hits_merged.bed',
                        help='BEDTools merge output')
    parser.add_argument('--output-tsv',
                        default='results/stage2/merged_regions.tsv',
                        help='Parsed merged regions with metadata')
    parser.add_argument('--merge-distance',
                        type=int, default=DEFAULT_MERGE_DISTANCE,
                        help='Merge hits within this many bp (default: 500)')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    print("\n==========================================")
    print(" Stage 2: Merge Nearby Hits")
    print("==========================================")

    # Step 1 — Check BEDTools
    check_bedtools()

    # Step 2 — Sort BED
    sort_bed(args.input_bed, args.sorted_bed)

    # Step 3 — Merge nearby hits
    merge_hits(args.sorted_bed, args.merged_bed, args.merge_distance)

    # Step 4 — Parse into clean TSV
    df = parse_merged_bed(args.merged_bed, args.output_tsv)

    # Step 5 — Summary
    print_summary(args.input_bed, df, args.merge_distance)

    print(f"\n==========================================")
    print(f" Stage 2 merge_nearby_hits — Complete")
    print(f"==========================================")
    print(f"  Output : {args.output_tsv}")
    print(f"  Next   : extend_regions.py")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()