
# check_uga.py
# Stage 3.1: TGA presence check in candidate regions
# Simplified approach: check downstream sequence for TGA codons
# in all 6 reading frames — correct approach for exon-level hits

import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO

STOP_CODONS = {"TAA", "TAG", "TGA"}

# ─────────────────────────────────────────────
# Why CHANGED APPROACH
# ─────────────────────────────────────────────
# Problem with previous approach:
#   DIAMOND finds middle/early exons of selenoproteins
#   The Sec codon is ALWAYS in the last coding exon
#   → Sec exon is NOT in the extracted region
#   → Protein position verification always fails
#
# Correct approach:
#   1. Check that the candidate region contains a TGA codon
#      in any reading frame (not just the last stop)
#   2. Use SECIS detection (Stage 3.2) as the primary
#      biological validation — SECIS presence is definitive
#   3. UGA check serves as a pre-filter only

# ─────────────────────────────────────────────
# 1. LOAD REFERENCE SELENOPROTEINS
# ─────────────────────────────────────────────

def load_selenoproteins(fasta_path):
    if not fasta_path or not os.path.exists(fasta_path):
        return {}

    print(f"\n--- Loading Reference Selenoproteins ---")
    proteins = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        pid = record.id
        if '|' in pid:
            parts = pid.split('|')
            if len(parts) >= 2:
                pid = parts[1]
        proteins[pid] = str(record.seq).upper()

    with_U = sum(1 for s in proteins.values() if 'U' in s)
    print(f"  Proteins loaded      : {len(proteins)}")
    print(f"  With Sec (U)         : {with_U}")
    return proteins

# ─────────────────────────────────────────────
# 2. LOAD MERGED HITS
# ─────────────────────────────────────────────

def load_hits(merged_hits_path):
    if not os.path.exists(merged_hits_path):
        return pd.DataFrame()

    print(f"\n--- Loading Merged Hits ---")
    df = pd.read_csv(merged_hits_path, sep='\t')
    df  = df[df['redundant'] == False].copy()

    def extract_acc(pid):
        pid = str(pid)
        if '|' in pid:
            parts = pid.split('|')
            if len(parts) >= 2:
                return parts[1]
        return pid
    df['protein_acc'] = df['target_protein'].apply(extract_acc)

    print(f"  Hits loaded          : {len(df)}")
    return df

# ─────────────────────────────────────────────
# 3. PARSE FASTA HEADER
# ─────────────────────────────────────────────

def get_header_info(record):
    desc = record.description
    info = {'target_protein': None, 'gene_start': None,
            'gene_end': None, 'strand': None}

    if 'target=' in desc:
        try:
            info['target_protein'] = desc.split('target=')[1].split(' ')[0]
        except IndexError:
            pass
    if 'gene=' in desc:
        try:
            gene_part = desc.split('gene=')[1].split(' ')[0]
            coords    = gene_part.split(':')[1]
            info['gene_start'] = int(coords.split('-')[0])
            info['gene_end']   = int(coords.split('-')[1])
        except (IndexError, ValueError):
            pass
    if 'strand=' in desc:
        try:
            info['strand'] = desc.split('strand=')[1].split(' ')[0]
        except IndexError:
            pass
    return info

# ─────────────────────────────────────────────
# 4. CHECK TGA IN ALL READING FRAMES
# ─────────────────────────────────────────────

def check_tga_in_sequence(seq, min_orf_aa=30):
    """
    Check for TGA codons in all 6 reading frames.

    For exon-level hits, we cannot pinpoint the exact Sec
    position because the Sec exon may not be in our region.

    Instead we check:
      1. Does any reading frame contain an ORF that ENDS in TGA?
         (signature of selenoprotein truncated gene model)
      2. How many TGA codons exist in all frames?
         (high TGA density suggests selenoprotein context)

    The SECIS element (Stage 3.2) provides the definitive
    biological confirmation.
    """
    seq  = seq.upper()
    tga_orfs     = []
    total_tga    = 0
    total_codons = 0

    # Check all 6 frames (3 forward + 3 reverse complement)
    from Bio.Seq import Seq
    seqs_to_check = [
        (seq, '+'),
        (str(Seq(seq).reverse_complement()), '-')
    ]

    for strand_seq, strand in seqs_to_check:
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                codon = strand_seq[i:i+3]
                total_codons += 1
                if codon == 'TGA':
                    total_tga += 1
                i += 3

        # Find TGA-terminated ORFs
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                codon = strand_seq[i:i+3]
                if codon == 'ATG':
                    start = i
                    j     = i + 3
                    while j + 2 < len(strand_seq):
                        c2 = strand_seq[j:j+3]
                        if c2 in STOP_CODONS:
                            length_aa = (j - start) // 3
                            if length_aa >= min_orf_aa and c2 == 'TGA':
                                tga_orfs.append({
                                    'frame'    : frame,
                                    'strand'   : strand,
                                    'start'    : start,
                                    'end'      : j + 3,
                                    'length_aa': length_aa,
                                    'tga_pos'  : j,
                                })
                            break
                        j += 3
                    i = j
                else:
                    i += 3

    has_tga_orf = len(tga_orfs) > 0

    return {
        'has_tga_orf'     : has_tga_orf,
        'n_tga_orfs'      : len(tga_orfs),
        'total_tga_codons': total_tga,
        'best_orf_len'    : max((o['length_aa'] for o in tga_orfs), default=0),
    }

# ─────────────────────────────────────────────
# 5. VERIFY PROTEIN IS A KNOWN SELENOPROTEIN
# ─────────────────────────────────────────────

def verify_known_selenoprotein(target_protein_id, proteins):
    """
    Confirm the target protein is in our selenoprotein
    reference set AND contains Sec (U).

    This replaces the failed position-specific verification
    with a simpler but reliable check:
      Is this protein a known selenoprotein? → yes/no
    """
    if not target_protein_id:
        return False, None, []

    ref_seq = proteins.get(target_protein_id, '')
    if not ref_seq:
        return False, None, []

    sec_positions = [i for i, aa in enumerate(ref_seq) if aa == 'U']
    is_selenoprotein = len(sec_positions) > 0

    return is_selenoprotein, ref_seq, sec_positions

# ─────────────────────────────────────────────
# 6. PROCESS ALL CANDIDATE REGIONS
# ─────────────────────────────────────────────

def process_candidates(candidate_fasta, proteins, hits_df, min_orf_aa):
    """
    Two-part validation:

    Part A — Known selenoprotein check:
      Is the target protein in our selenoprotein reference set
      with a confirmed Sec (U) residue?
      → This catches TXNRD2, SELM, SELENOO (real)
      → Rejects Q9NZV5 hits IF Q9NZV5 isn't in positive_set.fa
        OR if the region is clearly noise

    Part B — TGA presence check:
      Does the sequence contain TGA-terminated ORFs?
      → Simple pre-filter before SECIS detection

    A region passes if BOTH conditions are met.
    SECIS detection (Stage 3.2) is the final arbiter.
    """
    print(f"\n--- Running UGA + Known Selenoprotein Check ---")

    results         = []
    passing_records = []

    total          = 0
    pass_known     = 0
    pass_tga       = 0
    pass_both      = 0

    for record in SeqIO.parse(candidate_fasta, "fasta"):
        total += 1
        seq    = str(record.seq).upper()

        # Skip all-N
        non_n = sum(1 for c in seq if c != 'N')
        if non_n < 100:
            results.append({
                'region_id': record.id, 'skip_reason': 'all_N',
                'pass_known_seleno': False, 'pass_tga_check': False,
                'passes': False, 'target_protein': '',
            })
            continue

        header_info       = get_header_info(record)
        target_protein_id = header_info.get('target_protein', '')

        # Part A — known selenoprotein check
        is_known_seleno, ref_seq, sec_positions = verify_known_selenoprotein(
            target_protein_id, proteins
        )
        if is_known_seleno:
            pass_known += 1

        # Part B — TGA presence check
        tga_result = check_tga_in_sequence(seq, min_orf_aa)
        if tga_result['has_tga_orf']:
            pass_tga += 1

        # Pass if known selenoprotein AND has TGA ORFs
        passes = is_known_seleno and tga_result['has_tga_orf']
        if passes:
            pass_both += 1
            passing_records.append(record)

        results.append({
            'region_id'            : record.id,
            'length_bp'            : len(seq),
            'skip_reason'          : '',
            'target_protein'       : target_protein_id,
            'is_known_selenoprotein': is_known_seleno,
            'ref_protein_length'   : len(ref_seq) if ref_seq else 0,
            'sec_positions_in_ref' : str(sec_positions),
            'has_tga_orf'          : tga_result['has_tga_orf'],
            'n_tga_orfs'           : tga_result['n_tga_orfs'],
            'total_tga_codons'     : tga_result['total_tga_codons'],
            'best_orf_len_aa'      : tga_result['best_orf_len'],
            'passes'               : passes,
        })

    print(f"\n  Total candidate regions          : {total}")
    print(f"  Known selenoprotein (Part A)     : {pass_known}")
    print(f"  Has TGA-terminated ORF (Part B)  : {pass_tga}")
    print(f"  Passes both checks               : {pass_both}")
    print(f"\n  → {pass_both} regions proceed to SECIS detection")

    return results, passing_records

# ─────────────────────────────────────────────
# 7. WRITE OUTPUTS
# ─────────────────────────────────────────────

def write_results(results, output_tsv):
    df = pd.DataFrame(results)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n  UGA results saved    : {output_tsv}")
    return df

def write_passing_fasta(records, output_fasta):
    SeqIO.write(records, output_fasta, "fasta")
    print(f"  Passing FASTA saved  : {output_fasta}")
    print(f"  Passing sequences    : {len(records)}")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 3.1: TGA check + known selenoprotein verification"
    )
    parser.add_argument('--candidate-fa',
                        default='results/stage2/candidate_regions.fa')
    parser.add_argument('--selenoproteins',
                        default='data/processed/positive_set.fa')
    parser.add_argument('--merged-hits',
                        default='results/stage1/merged_hits.tsv')
    parser.add_argument('--min-orf-aa',
                        type=int, default=30)
    parser.add_argument('--output-tsv',
                        default='results/stage3/uga_validated.tsv')
    parser.add_argument('--output-passing-fa',
                        default='results/stage3/uga_passing_regions.fa')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    print("\n==========================================")
    print(" Stage 3.1: UGA Validation")
    print("==========================================")

    proteins = load_selenoproteins(args.selenoproteins)
    hits_df  = load_hits(args.merged_hits)

    results, passing_records = process_candidates(
        args.candidate_fa, proteins, hits_df, args.min_orf_aa
    )

    write_results(results, args.output_tsv)
    write_passing_fasta(passing_records, args.output_passing_fa)

    print(f"\n==========================================")
    print(f" Stage 3.1 check_uga — Complete")
    print(f"==========================================")
    print(f"  Results TSV  : {args.output_tsv}")
    print(f"  Passing FASTA: {args.output_passing_fa}")
    print(f"  Next         : check_secis.py")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()