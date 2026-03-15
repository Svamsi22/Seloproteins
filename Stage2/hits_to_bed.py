
# hits_to_bed.py
# Stage 2: Convert alignment hits to BED format
# KEY: Select most 3'-proximal hit per gene as extension anchor

import argparse
import os
import sys
import pandas as pd

def load_hits(merged_tsv):
    print(f"\n--- Loading Merged Hits ---")

    if not os.path.exists(merged_tsv):
        print(f"  ERROR: {merged_tsv} not found.")
        sys.exit(1)

    df = pd.read_csv(merged_tsv, sep='\t')

    total      = len(df)
    redundant  = df['redundant'].sum()
    to_process = df[df['redundant'] == False]

    print(f"  Total hits loaded    : {total}")
    print(f"  Redundant (skipped)  : {redundant}")
    print(f"  To process           : {len(to_process)}")

    return to_process

# ─────────────────────────────────────────────
# NORMALISE PROTEIN IDs
# ─────────────────────────────────────────────

def normalise_ids(df):
    def extract_acc(pid):
        if '|' in str(pid):
            return str(pid).split('|')[1]
        return str(pid)
    df = df.copy()
    df['protein_acc'] = df['target_protein'].apply(extract_acc)
    return df

# ─────────────────────────────────────────────
# KEY FIX: SELECT MOST 3'-PROXIMAL HIT PER GENE
# ─────────────────────────────────────────────

def select_anchor_hits(df):
    """
    For SECIS detection we MUST extend from the most 3'-proximal
    exon of each gene — not from random middle exons.

    Why:
      SECIS element is in the 3' UTR, after the LAST coding exon.
      If we extend from an early exon, we search the wrong region.

    Rule:
      + strand gene: most 3' exon = hit with HIGHEST query_end
      - strand gene: most 3' exon = hit with LOWEST query_start

    Example:
      SELM (- strand, 12 exon hits):
        Early exon hit:  31,107,376 - 31,107,505  (pos 1-43/145)
        Last  exon hit:  31,104,973 - 31,105,131  (pos 93-145/145) ← use this
        Extension from last exon → reaches SECIS in 3' UTR

    We group by protein_acc and keep one anchor per gene.
    """
    print(f"\n--- Selecting Most 3'-Proximal Exon Per Gene ---")

    selected = []

    for protein_acc, group in df.groupby('protein_acc'):
        strand = group['strand'].iloc[0]

        if strand == '+':
            anchor = group.loc[group['query_end'].idxmax()]
        else:
            anchor = group.loc[group['query_start'].idxmin()]

        selected.append(anchor)

        # Show selection for debugging
        print(f"  {protein_acc:<15} ({strand})  "
              f"anchor: {int(anchor['query_start'])}-{int(anchor['query_end'])}  "
              f"protein: {int(anchor['target_start'])}-{int(anchor['target_end'])}"
              f"/{int(anchor['target_length'])}  "
              f"bitscore: {anchor['bitscore']}")

    anchored_df = pd.DataFrame(selected).reset_index(drop=True)

    print(f"\n  Hits before anchor selection : {len(df)}")
    print(f"  Hits after  anchor selection : {len(anchored_df)}")

    return anchored_df

# ─────────────────────────────────────────────
# CONVERT TO BED
# ─────────────────────────────────────────────

def hits_to_bed(df):
    print(f"\n--- Converting Hits to BED Format ---")

    bed_rows = []

    for _, hit in df.iterrows():
        chrom  = str(hit['query_id'])
        start  = int(hit['query_start']) - 1   # 0-based
        end    = int(hit['query_end'])
        score  = int(hit['bitscore'])
        strand = hit['strand']

        protein_acc = hit.get('protein_acc', str(hit['target_protein']))
        name = (
            f"{hit['query_id']}|"
            f"{protein_acc}|"
            f"{hit['source']}|"
            f"{hit['identity']}|"
            f"{hit['coverage']}"
        )

        bed_rows.append({
            'chrom'  : chrom,
            'start'  : start,
            'end'    : end,
            'name'   : name,
            'score'  : score,
            'strand' : strand,
        })

    bed_df = pd.DataFrame(bed_rows)
    bed_df = bed_df.sort_values(['chrom', 'start']).reset_index(drop=True)

    print(f"  BED entries created  : {len(bed_df)}")
    print(f"  Chromosomes          : {bed_df['chrom'].unique().tolist()}")
    return bed_df

# ─────────────────────────────────────────────
# VALIDATE CHROMOSOME NAMES
# ─────────────────────────────────────────────

def validate_chroms(bed_df, genome_fai):
    print(f"\n--- Validating Chromosome Names ---")

    if not os.path.exists(genome_fai):
        print(f"  WARNING: .fai not found — skipping check")
        return

    fai_chroms = set()
    with open(genome_fai) as f:
        for line in f:
            fai_chroms.add(line.split('\t')[0])

    bed_chroms     = set(bed_df['chrom'].unique())
    missing_chroms = bed_chroms - fai_chroms

    if missing_chroms:
        print(f"  ❌ Chromosome mismatch: {missing_chroms}")
        sys.exit(1)
    else:
        print(f" All chromosome names match genome index")

# ─────────────────────────────────────────────
# WRITE BED
# ─────────────────────────────────────────────

def write_bed(bed_df, output_bed):
    bed_df.to_csv(output_bed, sep='\t', index=False, header=False)
    print(f"\n  BED file saved       : {output_bed}")
    print(f"  Entries              : {len(bed_df)}")

    print(f"\n  All entries:")
    for _, row in bed_df.iterrows():
        print(f"    {row['chrom']}  {row['start']}  {row['end']}  "
              f"{row['name'].split('|')[1]:<15}  "
              f"score={row['score']}  strand={row['strand']}")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: Convert alignment hits to BED "
                    "(anchored at most 3'-proximal exon per gene)"
    )
    parser.add_argument('--merged-hits',
                        default='results/stage1/merged_hits.tsv')
    parser.add_argument('--output-bed',
                        default='results/stage2/hits.bed')
    parser.add_argument('--genome-fai',
                        default='data/raw/chr22_clean.fa.fai')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_bed), exist_ok=True)

    print("\n==========================================")
    print(" Stage 2: Hits to BED (3' anchor selection)")
    print("==========================================")

    # Step 1 — Load non-redundant hits
    df = load_hits(args.merged_hits)

    # Step 2 — Normalise protein IDs
    df = normalise_ids(df)

    # Step 3 — Select most 3'-proximal hit per gene
    df = select_anchor_hits(df)

    # Step 4 — Convert to BED
    bed_df = hits_to_bed(df)

    # Step 5 — Validate chromosome names
    validate_chroms(bed_df, args.genome_fai)

    # Step 6 — Write BED
    write_bed(bed_df, args.output_bed)

    print(f"\n==========================================")
    print(f" Stage 2 hits_to_bed — Complete")
    print(f"==========================================")
    print(f"  Output : {args.output_bed}")
    print(f"  Next   : merge_nearby_hits.py")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()