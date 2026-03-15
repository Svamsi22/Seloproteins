# extract_ground_truth.py
# Extracts selenoprotein entries from full Ensembl GFF3
# Uses multiple detection strategies to ensure nothing is missed

import argparse
import re
import sys
from collections import defaultdict

# ─────────────────────────────────────────────
# KNOWN SELENOPROTEINS ON CHROMOSOME 22
# ─────────────────────────────────────────────
# These are experimentally verified selenoproteins
# known to be encoded on human chromosome 22
# Source: NCBI Gene + Ensembl annotation

KNOWN_SELENO_GENE_NAMES = {
    'TXNRD2',    # Thioredoxin reductase 2 (mitochondrial)
    'SELENOO',   # Selenoprotein O
    'GPX6',      # Glutathione peroxidase 6 (pseudogene in some individuals)
    'SELENOM',   # Selenoprotein M (check if on chr22)
}

# ─────────────────────────────────────────────
# 1. PARSE FULL GFF3 INTO MEMORY
# ─────────────────────────────────────────────

def parse_gff3(gff3_path):
    """
    Parse GFF3 into:
      - header lines (comments)
      - features dict keyed by feature ID
      - parent->children mapping
    """
    print(f"\n  Reading: {gff3_path}")

    headers      = []
    features     = {}       # ID -> full line
    children     = defaultdict(list)   # parent_id -> [child_ids]
    feature_list = []       # preserve order

    with open(gff3_path) as f:
        for line in f:

            # Keep header lines
            if line.startswith('#'):
                headers.append(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            attributes = parts[8]

            # Extract ID
            id_match = re.search(r'(?:^|;)ID=([^;]+)', attributes)
            feat_id  = id_match.group(1) if id_match else None

            # Extract Parent
            parent_match = re.search(r'(?:^|;)Parent=([^;]+)', attributes)
            parent_ids   = parent_match.group(1).split(',') \
                           if parent_match else []

            if feat_id:
                features[feat_id] = line
                feature_list.append(feat_id)

            # Build parent → children map
            for pid in parent_ids:
                if feat_id:
                    children[pid].append(feat_id)

    total = len(features)
    print(f"  Total features parsed : {total:,}")
    return headers, features, children, feature_list

# ─────────────────────────────────────────────
# 2. DETECT SELENOPROTEIN GENES
# ─────────────────────────────────────────────

def detect_selenoprotein_genes(features):
    """
    Two strategies to find selenoprotein gene IDs:

    Strategy 1: Gene name matches known selenoprotein list
    Strategy 2: Any feature has 'selenocysteine' in attributes
    """
    found_by_name  = set()
    found_by_annot = set()

    for feat_id, line in features.items():
        parts      = line.strip().split('\t')
        if len(parts) < 9:
            continue

        feat_type  = parts[2]
        attributes = parts[8]

        # ── Strategy 1: Known gene names ──────────────────────
        if feat_type == 'gene':
            name_match = re.search(r'(?:^|;)Name=([^;]+)', attributes)
            if name_match:
                gene_name = name_match.group(1).upper()
                if gene_name in KNOWN_SELENO_GENE_NAMES:
                    found_by_name.add(feat_id)
                    print(f"    [name match]  {feat_id} → {gene_name}")

        # ── Strategy 2: Selenocysteine annotation ─────────────
        if 'selenocysteine' in attributes.lower():
            # Walk up to find the parent gene ID
            parent_match = re.search(r'(?:^|;)Parent=([^;]+)', attributes)
            if parent_match:
                for pid in parent_match.group(1).split(','):
                    found_by_annot.add(pid)
                    print(f"    [annot match] {feat_id} → parent {pid}")
            elif feat_type == 'gene':
                found_by_annot.add(feat_id)
                print(f"    [annot match] {feat_id} (gene level)")

    return found_by_name, found_by_annot

# ─────────────────────────────────────────────
# 3. COLLECT ALL CHILD FEATURES
# ─────────────────────────────────────────────

def collect_all_children(gene_ids, children):
    """
    Given a set of gene IDs, recursively collect all
    descendant feature IDs (mRNA → CDS, exon, UTR etc.)

    GFF3 hierarchy:
      gene
       └── mRNA
            ├── CDS
            ├── exon
            └── UTR
    """
    all_ids = set(gene_ids)
    queue   = list(gene_ids)

    while queue:
        current = queue.pop()
        for child_id in children.get(current, []):
            if child_id not in all_ids:
                all_ids.add(child_id)
                queue.append(child_id)

    return all_ids

# ─────────────────────────────────────────────
# 4. WRITE OUTPUT GFF3
# ─────────────────────────────────────────────

def write_ground_truth_gff3(headers, features, feature_list,
                              keep_ids, output_path):
    """
    Write filtered GFF3 preserving original feature order.
    Only includes features whose ID is in keep_ids.
    """
    kept = 0
    with open(output_path, 'w') as out:

        # Write headers
        out.writelines(headers)
        out.write("## Filtered by extract_ground_truth.py\n")
        out.write("## Contains selenoprotein genes only\n")

        # Write features in original order
        for feat_id in feature_list:
            if feat_id in keep_ids:
                out.write(features[feat_id])
                kept += 1

    return kept

# ─────────────────────────────────────────────
# 5. PRINT SUMMARY
# ─────────────────────────────────────────────

def print_summary(found_by_name, found_by_annot,
                  all_gene_ids, total_features_kept):
    """
    Clear summary of what was found and how.
    """
    print("\n==========================================")
    print(" Ground Truth Extraction Summary")
    print("==========================================")
    print(f"  Found by gene name   : {len(found_by_name)}")
    print(f"  Found by annotation  : {len(found_by_annot)}")
    print(f"  Total unique genes   : {len(all_gene_ids)}")
    print(f"  Total features kept  : {total_features_kept}")

    if not all_gene_ids:
        print("\n  WARNING: No selenoproteins found in GFF3")
        print("  Possible reasons:")
        print("    - Chromosome naming mismatch")
        print("    - Gene names differ from expected")
        print("    - Wrong Ensembl release")
    else:
        print("\n  Selenoprotein genes extracted:")
        for gid in sorted(all_gene_ids):
            print(f"    {gid}")

    print("==========================================\n")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Extract selenoprotein ground truth from Ensembl GFF3"
    )
    parser.add_argument('--input',  default='data/ground_truth/ensembl_chr22_full.gff3',
                        help='Full Ensembl GFF3 (ensembl_chr22_full.gff3)')
    parser.add_argument('--output', default='data/ground_truth/ground_truth_selenoproteins.gff3',
                        help='Output filtered GFF3')
    args = parser.parse_args()

    print("\n==========================================")
    print(" extract_ground_truth.py")
    print("==========================================")

    # Step 1 — Parse full GFF3
    headers, features, children, feature_list = parse_gff3(args.input)

    # Step 2 — Detect selenoprotein genes
    print("\n  Detecting selenoprotein genes...")
    found_by_name, found_by_annot = detect_selenoprotein_genes(features)

    # Step 3 — Merge both detection sets
    all_gene_ids = found_by_name | found_by_annot

    if not all_gene_ids:
        print("\n  WARNING: No selenoprotein genes detected")
        print("  Writing empty GFF3 — check input file")

    # Step 4 — Collect all child features
    all_feature_ids = collect_all_children(all_gene_ids, children)
    print(f"\n  Total features to extract : {len(all_feature_ids)}")

    # Step 5 — Write filtered GFF3
    kept = write_ground_truth_gff3(
        headers, features, feature_list,
        all_feature_ids, args.output
    )

    # Step 6 — Summary
    print_summary(found_by_name, found_by_annot, all_gene_ids, kept)
    print(f"  Saved to: {args.output}")

if __name__ == '__main__':
    main()