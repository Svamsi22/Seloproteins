
# extend_regions.py
# Stage 2: Extend merged regions to include downstream sequence
# Downstream extension captures SECIS elements for Stage 3 validation

import argparse
import os
import sys
import pandas as pd

# ─────────────────────────────────────────────
# WHY WE EXTEND
# ─────────────────────────────────────────────
# SECIS elements are always in the 3' UTR — AFTER the stop codon.
# They can be up to 5kb downstream of the UGA codon in eukaryotes.
#
# Without extension:
#   Region: |===gene===|
#   SECIS:              |--secis--|   ← missed, outside region
#
# With extension:
#   Region: |===gene===|----buffer----|
#   SECIS:              |--secis--|   ← captured 
#
# Extension is STRAND-AWARE:
#   + strand gene: extend the END coordinate
#   - strand gene: extend the START coordinate
#
#   + strand:  |===gene===>----5kb----|
#   - strand:  |----5kb----<===gene===|

# ─────────────────────────────────────────────
# 1. LOAD GENOME CHROMOSOME SIZES FROM .FAI
# ─────────────────────────────────────────────

def load_chrom_sizes(genome_fai):
    """
    Read chromosome sizes from samtools .fai index.
    We need these to make sure extended regions
    don't go beyond the end of the chromosome.

    .fai format:
      chrom   length   offset   bases_per_line   bytes_per_line
    """
    print(f"\n--- Loading Chromosome Sizes ---")

    if not os.path.exists(genome_fai):
        print(f"  ERROR: .fai not found: {genome_fai}")
        print(f"  Run: samtools faidx data/raw/chr22_clean.fa")
        sys.exit(1)

    chrom_sizes = {}
    with open(genome_fai) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])

    for chrom, size in chrom_sizes.items():
        print(f"  {chrom:<6} : {size:,} bp  ({size/1e6:.1f} Mb)")

    return chrom_sizes

# ─────────────────────────────────────────────
# 2. EXTEND A SINGLE REGION
# ─────────────────────────────────────────────

def extend_region(row, chrom_sizes, upstream_buffer, extension_buffer):
    """
    Extend one region based on strand.

    upstream_buffer  : small extension at gene start
                       captures any missed start codons (default 200bp)

    extension_buffer : large extension at gene end
                       captures SECIS element (default 5000bp)

    Strand logic:
      + strand: gene runs left → right
                downstream = higher coordinates
                extend END by extension_buffer
                extend START back by upstream_buffer

      - strand: gene runs right → left
                downstream = lower coordinates
                extend START by extension_buffer
                extend END forward by upstream_buffer

    Boundary clamping:
      Never go below 0 or above chromosome length.
    """
    chrom    = row['chrom']
    start    = row['start']
    end      = row['end']
    strand   = row['strand']
    chrom_len = chrom_sizes.get(str(chrom), 0)

    if chrom_len == 0:
        print(f"  WARNING: chrom {chrom} not found in .fai — skipping")
        return None

    if strand == '+':
        ext_start = max(0, start - upstream_buffer)
        ext_end   = min(chrom_len, end + extension_buffer)
    else:
        ext_start = max(0, start - extension_buffer)
        ext_end   = min(chrom_len, end + upstream_buffer)

    # Track how much was actually added
    # (may be less than buffer if near chromosome boundary)
    actual_upstream = start - ext_start
    actual_downstream = ext_end - end

    return {
        'ext_start'         : ext_start,
        'ext_end'           : ext_end,
        'actual_upstream'   : actual_upstream,
        'actual_downstream' : actual_downstream,
        'region_size_bp'    : ext_end - ext_start,
        'clamped'           : (
            ext_start == 0 or ext_end == chrom_len
        ),
    }

# ─────────────────────────────────────────────
# 3. EXTEND ALL REGIONS
# ─────────────────────────────────────────────

def extend_all_regions(df, chrom_sizes, upstream_buffer, extension_buffer):
    """
    Apply extension to every merged region.
    Adds extended coordinate columns alongside original ones.
    """
    print(f"\n--- Extending Regions ---")
    print(f"  Upstream buffer    : {upstream_buffer:,} bp")
    print(f"  Downstream buffer  : {extension_buffer:,} bp")
    print(f"  Regions to extend  : {len(df)}")

    extended_rows = []
    skipped       = 0
    clamped       = 0

    for _, row in df.iterrows():
        result = extend_region(
            row, chrom_sizes,
            upstream_buffer,
            extension_buffer
        )

        if result is None:
            skipped += 1
            continue

        if result['clamped']:
            clamped += 1
            print(f" Clamped at chromosome boundary: {row['region_id']}")

        extended_rows.append({
            **row.to_dict(),

            # Extended coordinates
            'ext_start'         : result['ext_start'],
            'ext_end'           : result['ext_end'],

            # Keep original gene coordinates for UGA checking in Stage 3
            'gene_start'        : row['start'],
            'gene_end'          : row['end'],

            # Extension metadata
            'actual_upstream'   : result['actual_upstream'],
            'actual_downstream' : result['actual_downstream'],
            'ext_region_size_bp': result['region_size_bp'],
            'clamped'           : result['clamped'],

            # Updated region ID using extended coordinates
            'ext_region_id'     : (
                f"{row['chrom']}:"
                f"{result['ext_start']}-"
                f"{result['ext_end']}"
                f"({row['strand']})"
            ),
        })

    ext_df = pd.DataFrame(extended_rows)

    print(f"\n  Regions extended     : {len(ext_df)}")
    print(f"  Skipped (no chrom)   : {skipped}")
    print(f"  Clamped (at boundary): {clamped}")

    return ext_df

# ─────────────────────────────────────────────
# 4. WRITE EXTENDED BED
# ─────────────────────────────────────────────

def write_extended_bed(ext_df, output_bed):
    """
    Write extended regions as BED file for bedtools getfasta.

    Uses EXTENDED coordinates (ext_start, ext_end).
    Name field encodes ext_region_id for tracking.

    BED format reminder:
      col 0: chrom
      col 1: ext_start   (0-based)
      col 2: ext_end
      col 3: ext_region_id
      col 4: max_bitscore
      col 5: strand
    """
    bed_rows = ext_df[[
        'chrom',
        'ext_start',
        'ext_end',
        'ext_region_id',
        'max_bitscore',
        'strand',
    ]].copy()

    bed_rows.to_csv(output_bed, sep='\t', index=False, header=False)
    print(f"\n  Extended BED saved   : {output_bed}")
    print(f"  Entries              : {len(bed_rows)}")

# ─────────────────────────────────────────────
# 5. SUMMARY
# ─────────────────────────────────────────────

def print_summary(df, ext_df):
    """
    Show before/after region sizes to confirm extension worked.
    """
    print(f"\n--- Extension Summary ---")
    print(f"  {'Region':<35} {'Before':>10}  {'After':>10}  "
          f"{'Strand':>7}  {'Clamped':>8}")
    print(f"  {'-'*35} {'-'*10}  {'-'*10}  {'-'*7}  {'-'*8}")

    for _, row in ext_df.iterrows():
        before = row['region_size_bp']
        after  = row['ext_region_size_bp']
        print(
            f"  {row['region_id']:<35} "
            f"{before:>10,}  "
            f"{after:>10,}  "
            f"{row['strand']:>7}  "
            f"{'yes' if row['clamped'] else 'no':>8}"
        )

    if ext_df.empty:
        return

    print(f"\n  Size after extension:")
    print(f"    Min  : {ext_df['ext_region_size_bp'].min():,} bp")
    print(f"    Max  : {ext_df['ext_region_size_bp'].max():,} bp")
    print(f"    Mean : {ext_df['ext_region_size_bp'].mean():,.0f} bp")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: Extend merged regions for SECIS capture"
    )
    parser.add_argument('--merged-regions',
                        default='results/stage2/merged_regions.tsv',
                        help='Merged regions from merge_nearby_hits.py')
    parser.add_argument('--genome-fai',
                        default='data/raw/chr22_clean.fa.fai',
                        help='Genome .fai index for chromosome sizes')
    parser.add_argument('--output-tsv',
                        default='results/stage2/extended_regions.tsv',
                        help='Extended regions with all metadata')
    parser.add_argument('--output-bed',
                        default='results/stage2/extended_regions.bed',
                        help='Extended regions BED for bedtools getfasta')
    parser.add_argument('--extension-buffer',
                        type=int, default=8000,
                        help='Downstream extension in bp (default: 8000)')
    parser.add_argument('--upstream-buffer',
                        type=int, default=300,
                        help='Upstream extension in bp (default: 300)')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    print("\n==========================================")
    print(" Stage 2: Extend Regions")
    print("==========================================")

    # Step 1 — Load chromosome sizes
    chrom_sizes = load_chrom_sizes(args.genome_fai)

    # Step 2 — Load merged regions
    print(f"\n--- Loading Merged Regions ---")
    df = pd.read_csv(args.merged_regions, sep='\t')
    print(f"  Regions loaded : {len(df)}")

    # Step 3 — Extend all regions
    ext_df = extend_all_regions(
        df, chrom_sizes,
        args.upstream_buffer,
        args.extension_buffer
    )

    if ext_df.empty:
        print("  ERROR: No regions after extension.")
        sys.exit(1)

    # Step 4 — Write outputs
    ext_df.to_csv(args.output_tsv, sep='\t', index=False)
    print(f"\n  Extended TSV saved   : {args.output_tsv}")

    write_extended_bed(ext_df, args.output_bed)

    # Step 5 — Summary
    print_summary(df, ext_df)

    print(f"\n==========================================")
    print(f" Stage 2 extend_regions — Complete")
    print(f"==========================================")
    print(f"  TSV output : {args.output_tsv}")
    print(f"  BED output : {args.output_bed}")
    print(f"  Next       : extract_sequences.py")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()