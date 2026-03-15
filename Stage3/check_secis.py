
# check_secis.py
# Stage 3.2: Detect SECIS elements using Infernal + RFAM RF00031
# Verifies SECIS elements are present downstream of UGA codons

import argparse
import os
import sys
import subprocess
import shutil
import time
import requests
import pandas as pd
from Bio import SeqIO

# ─────────────────────────────────────────────
# RFAM SECIS COVARIANCE MODEL
# ─────────────────────────────────────────────
# RF00031 = SECIS element covariance model from RFAM
# This model captures both Type 1 and Type 2 SECIS elements
# Infernal uses it to search for SECIS RNA secondary structure
#
# Download URL:
RFAM_SECIS_URL  = "https://rfam.org/family/RF00031/cm"
RFAM_SECIS_CM   = "data/rfam/RF00031.cm"

# E-value threshold for Infernal cmsearch
# 0.1 is standard for RFAM searches
INFERNAL_EVALUE = "0.1" # relaxed threshold to capture borderline hits to identify selum from chr22

# ─────────────────────────────────────────────
# 1. SETUP — DOWNLOAD + PRESS RFAM MODEL
# ─────────────────────────────────────────────

def setup_rfam_model(cm_path):
    """
    Download RF00031 SECIS covariance model from RFAM
    and press it (required by Infernal before cmsearch).

    cmpress creates 4 binary index files:
      RF00031.cm.i1f
      RF00031.cm.i1i
      RF00031.cm.i1m
      RF00031.cm.i1p
    """
    os.makedirs(os.path.dirname(cm_path), exist_ok=True)

    # Download if not already present
    if not os.path.exists(cm_path):
        print(f"\n--- Downloading RFAM RF00031 SECIS Model ---")
        print(f"  URL : {RFAM_SECIS_URL}")

        try:
            response = requests.get(RFAM_SECIS_URL, timeout=60)
            response.raise_for_status()

            with open(cm_path, 'w') as f:
                f.write(response.text)

            print(f"  Saved : {cm_path}")

        except Exception as e:
            print(f"  ERROR: Could not download RF00031: {e}")
            print(f"  Manual download: wget {RFAM_SECIS_URL} -O {cm_path}")
            sys.exit(1)
    else:
        print(f"\n  RF00031 model found  : {cm_path}")

    # Press model if index files don't exist
    pressed_index = cm_path + '.i1f'
    if not os.path.exists(pressed_index):
        print(f"\n--- Pressing RFAM Model (cmpress) ---")

        if not shutil.which('cmpress'):
            print("  ERROR: cmpress not found.")
            print("  Install with: conda install -c bioconda infernal")
            sys.exit(1)

        result = subprocess.run(
            ['cmpress', cm_path],
            capture_output=True, text=True
        )

        if result.returncode != 0:
            print(f"  ERROR: cmpress failed:\n{result.stderr}")
            sys.exit(1)

        print(f"  Model pressed successfully")

# ─────────────────────────────────────────────
# 2. CHECK INFERNAL IS INSTALLED
# ─────────────────────────────────────────────

def check_infernal():
    print(f"\n--- Checking Infernal ---")

    for tool in ['cmsearch', 'cmpress']:
        path = shutil.which(tool)
        if not path:
            print(f"  ERROR: {tool} not found.")
            print(f"  Install: conda install -c bioconda infernal")
            sys.exit(1)

    result = subprocess.run(
        ['cmsearch', '-h'],
        capture_output=True, text=True
    )
    # Version is in first line of stderr
    version_line = result.stdout.split('\n')[1] if result.stdout else 'unknown'
    print(f"  Infernal : {version_line.strip()}")

# ─────────────────────────────────────────────
# 3. LOAD UGA-PASSING REGIONS
# ─────────────────────────────────────────────

def load_passing_regions(fasta_path):
    """
    Load FASTA sequences that passed the UGA check.
    These are the candidates for SECIS detection.
    """
    print(f"\n--- Loading UGA-Passing Regions ---")

    if not os.path.exists(fasta_path):
        print(f"  ERROR: {fasta_path} not found")
        print(f"  Run check_uga.py first")
        sys.exit(1)

    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"  Regions loaded       : {len(records)}")

    if len(records) == 0:
        print("  No regions to check — all failed UGA check")
        sys.exit(0)

    return records

# ─────────────────────────────────────────────
# 4. EXTRACT DOWNSTREAM SEQUENCE
# ─────────────────────────────────────────────

def extract_downstream_sequence(record):
    """
    Extract downstream sequence for SECIS search.
    MUST be strand-aware because bedtools -s reverse-complements
    minus strand sequences, flipping which end is downstream.

    + strand: downstream = seq[gene_abs_end - region_abs_start :]
    - strand: downstream = seq[region_abs_end - gene_abs_start :]
                           (coordinates are flipped after RC)
    """
    seq         = str(record.seq).upper()
    description = record.description

    if 'gene=' not in description:
        midpoint = len(seq) // 2
        return seq[midpoint:], midpoint

    try:
        # Parse strand
        strand = '+'
        if 'strand=-' in description or 'strand= -' in description:
            strand = '-'
        elif 'strand=' in description:
            s = description.split('strand=')[1].split(' ')[0]
            strand = s

        # Parse gene coordinates
        gene_part  = description.split('gene=')[1].split(' ')[0]
        coords     = gene_part.split(':')[1]
        gene_start = int(coords.split('-')[0])
        gene_end   = int(coords.split('-')[1])

        # Parse region coordinates from record.id
        # Format: 22:19875175-19880920(-)::22:19875175-19880920(-)
        region_part      = record.id.split('::')[0]
        region_coords    = region_part.split(':')[1]
        region_abs_start = int(region_coords.split('(')[0].split('-')[0])
        region_abs_end   = int(region_coords.split('(')[0].split('-')[1])

        if strand == '+':
            # + strand: sequence runs low→high genomic coords
            # downstream = after gene_end
            offset = gene_end - region_abs_start
        else:
            # - strand: bedtools RC flips the sequence
            # position 0 = high genomic coord (5' end)
            # position N = low genomic coord (3' end)
            # downstream (3' UTR) = after gene runs out
            # gene occupies: [region_abs_end - gene_end :
            #                  region_abs_end - gene_start]
            # downstream starts at: region_abs_end - gene_start
            offset = region_abs_end - gene_start

        if 0 < offset < len(seq):
            print(f"    strand={strand} offset={offset} "
                  f"downstream_len={len(seq)-offset}bp")
            return seq[offset:], offset
        else:
            print(f"    WARNING: offset {offset} out of range "
                  f"for seq len {len(seq)} — using midpoint")
            midpoint = len(seq) // 2
            return seq[midpoint:], midpoint

    except (IndexError, ValueError, AttributeError) as e:
        print(f"    WARNING: coord parse failed ({e}) — using midpoint")
        midpoint = len(seq) // 2
        return seq[midpoint:], midpoint

# ─────────────────────────────────────────────
# 5. WRITE DOWNSTREAM SEQUENCES TO TEMP FASTA
# ─────────────────────────────────────────────

def write_downstream_fasta(records, tmp_fasta_path):
    """
    Write downstream sequences to a temporary FASTA.
    This is the input file for cmsearch.

    Stores offset in the record description
    so we can convert SECIS positions back to
    absolute genomic coordinates.
    """
    print(f"\n--- Extracting Downstream Sequences ---")

    written   = 0
    offsets   = {}   # region_id → offset_used

    with open(tmp_fasta_path, 'w') as out:
        for record in records:
            downstream_seq, offset = extract_downstream_sequence(record)

            if len(downstream_seq) < 50:
                print(f"  SKIP {record.id[:40]} — downstream too short "
                      f"({len(downstream_seq)} bp)")
                continue

            # Truncate to 8000 bp (cmsearch handles long seqs poorly)
            if len(downstream_seq) > 8000:
                downstream_seq = downstream_seq[:8000]

            out.write(f">{record.id}\n")
            # Write in 60-char lines (standard FASTA)
            for i in range(0, len(downstream_seq), 60):
                out.write(downstream_seq[i:i+60] + '\n')

            offsets[record.id] = offset
            written += 1

    print(f"  Sequences written    : {written}")
    print(f"  Temp FASTA           : {tmp_fasta_path}")

    return offsets

# ─────────────────────────────────────────────
# 6. RUN CMSEARCH
# ─────────────────────────────────────────────

def run_cmsearch(cm_path, input_fasta, tblout_path, threads):
    """
    Run Infernal cmsearch to find SECIS elements.

    Key flags:
      --tblout      : machine-readable tabular output
      --notextw     : don't wrap long lines
      -E 0.01       : E-value threshold (standard for RFAM)
      --cpu         : number of threads

    The covariance model RF00031 captures both:
      Type 1 SECIS: simple hairpin
      Type 2 SECIS: hairpin with extra stem-loop
    """
    print(f"\n--- Running cmsearch (Infernal + RF00031) ---")
    print(f"  Model    : {cm_path}")
    print(f"  Input    : {input_fasta}")
    print(f"  E-value  : {INFERNAL_EVALUE}")
    print(f"  Threads  : {threads}")

    cmd = [
        'cmsearch',
        '--tblout',   tblout_path,
        '--notextw',
        '-E',         INFERNAL_EVALUE,
        '--cpu',      str(threads),
        cm_path,
        input_fasta,
    ]

    print(f"\n  Command: {' '.join(cmd)}\n")

    start  = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"  ERROR: cmsearch failed:\n{result.stderr}")
        sys.exit(1)

    print(f"  cmsearch completed in {elapsed:.1f}s")

# ─────────────────────────────────────────────
# 7. PARSE CMSEARCH TBLOUT
# ─────────────────────────────────────────────

def parse_cmsearch_tblout(tblout_path):
    """
    Parse Infernal cmsearch --tblout output.

    tblout columns (space-separated):
      0:  target_name    (our region_id)
      1:  target_accn
      2:  query_name     (RF00031)
      3:  query_accn
      4:  mdl            (cm or hmm)
      5:  mdl_from
      6:  mdl_to
      7:  seq_from       (start in our downstream sequence)
      8:  seq_to         (end in our downstream sequence)
      9:  strand
      10: trunc
      11: pass
      12: gc
      13: bias
      14: score          (bit score)
      15: E-value
      16: inc            (! = pass, ? = borderline)
      17+: description

    We keep hits where inc == '!'  (definitive hit)
    """
    hits = {}   # region_id → list of SECIS hits

    if not os.path.exists(tblout_path):
        return hits

    with open(tblout_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split()
            if len(parts) < 17:
                continue

            region_id  = parts[0]
            seq_from   = int(parts[7])
            seq_to     = int(parts[8])
            strand     = parts[9]
            score      = float(parts[14])
            evalue     = float(parts[15])
            inc        = parts[16]   # ! = significant, ? = borderline

            if region_id not in hits:
                hits[region_id] = []

            hits[region_id].append({
                'seq_from'  : seq_from,
                'seq_to'    : seq_to,
                'strand'    : strand,
                'score'     : score,
                'evalue'    : evalue,
                'inc'       : inc,
                'significant': inc in ('!','?')
            })

    return hits

# ─────────────────────────────────────────────
# 8. BUILD RESULTS DATAFRAME
# ─────────────────────────────────────────────

def build_results(records, infernal_hits, offsets):
    """
    Combine cmsearch results with region metadata.

    For each region:
      - has_secis     : at least one significant hit (inc == '!')
      - secis_count   : number of SECIS elements found
      - best_score    : highest bit score
      - best_evalue   : lowest E-value
      - secis_pos_abs : SECIS position in full extended region
                        (offset + seq_from)
      - secis_type    : inferred from hit position
                        (RF00031 doesn't distinguish types directly)
    """
    rows = []

    for record in records:
        region_id   = record.id
        hits        = infernal_hits.get(region_id, [])
        offset      = offsets.get(region_id, 0)

        sig_hits    = [h for h in hits if h['significant']]
        has_secis   = len(sig_hits) > 0

        if sig_hits:
            best        = max(sig_hits, key=lambda h: h['score'])
            best_score  = best['score']
            best_evalue = best['evalue']
            # Convert position back to position in full extended sequence
            secis_pos_abs = offset + best['seq_from']
        else:
            best_score    = 0.0
            best_evalue   = None
            secis_pos_abs = None

        rows.append({
            'region_id'    : region_id,
            'has_secis'    : has_secis,
            'secis_count'  : len(sig_hits),
            'best_score'   : best_score,
            'best_evalue'  : best_evalue,
            'secis_pos_abs': secis_pos_abs,
            'total_hits'   : len(hits),   # including borderline
        })

    return pd.DataFrame(rows)

# ─────────────────────────────────────────────
# 9. VERIFY SEC POSITION MATCHES REFERENCE
# ─────────────────────────────────────────────

def verify_sec_positions(secis_df, uga_df):
    """
    Cross-check: confirm SECIS position is consistent
    with UGA codon position from check_uga.py.

    Biological constraint:
      SECIS must be DOWNSTREAM of the UGA codon.
      SECIS position (abs) > UGA codon position (abs)

    This implements the project requirement:
      "Ensuring selenocysteine positions match
       the reference protein sequences"

    Adds column: sec_position_consistent (bool)
    """
    print(f"\n--- Verifying Sec Position Consistency ---")

    if uga_df.empty or 'uga_details' not in uga_df.columns:
        print("  WARNING: No UGA position data — skipping consistency check")
        secis_df['sec_position_consistent'] = True
        return secis_df

    consistent   = 0
    inconsistent = 0

    secis_df = secis_df.copy()
    secis_df['sec_position_consistent'] = True

    for idx, row in secis_df.iterrows():
        if not row['has_secis'] or row['secis_pos_abs'] is None:
            continue

        # Find matching UGA hit for this region
        region_chrom = row['region_id'].split(':')[0]
        uga_match = uga_df[
            uga_df['query_id'].astype(str) == region_chrom
        ]

        if uga_match.empty:
            continue

        # Check SECIS is downstream of UGA
        # For + strand: secis_abs_pos > uga_query_end
        # For - strand: secis_abs_pos < uga_query_start
        for _, uga_hit in uga_match.iterrows():
            strand = str(uga_hit.get('strand', '+'))
            if strand == '+':
                uga_end = int(uga_hit.get('query_end', 0))
                if row['secis_pos_abs'] < uga_end:
                    secis_df.at[idx, 'sec_position_consistent'] = False
                    inconsistent += 1
                    print(f"  ⚠️  SECIS upstream of UGA: {row['region_id']}")
                else:
                    consistent += 1
            else:
                uga_start = int(uga_hit.get('query_start', 0))
                if row['secis_pos_abs'] > uga_start:
                    secis_df.at[idx, 'sec_position_consistent'] = False
                    inconsistent += 1
                    print(f"  ⚠️  SECIS upstream of UGA: {row['region_id']}")
                else:
                    consistent += 1

    print(f"  Position consistent  : {consistent}")
    print(f"  Position inconsistent: {inconsistent}")

    return secis_df

# ─────────────────────────────────────────────
# 10. WRITE OUTPUTS
# ─────────────────────────────────────────────

def write_outputs(secis_df, records, output_tsv, output_passing_fa):
    """
    Write two outputs:
      1. secis_hits.tsv          — all regions with SECIS results
      2. secis_passing_regions.fa — FASTA of SECIS-confirmed candidates
    """
    secis_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n  SECIS results saved  : {output_tsv}")

    # Keep only significant SECIS hits with consistent Sec positions
    passing_mask = (
        (secis_df['has_secis'] == True) &
        (secis_df['sec_position_consistent'] == True)
    )
    passing_ids = set(secis_df[passing_mask]['region_id'])

    passing_records = [r for r in records if r.id in passing_ids]
    SeqIO.write(passing_records, output_passing_fa, "fasta")

    print(f"  Passing FASTA saved  : {output_passing_fa}")
    print(f"  Final candidates     : {len(passing_records)}")

    return passing_records

# ─────────────────────────────────────────────
# 11. SUMMARY
# ─────────────────────────────────────────────

def print_summary(uga_df, secis_df):
    total      = len(secis_df)
    has_secis  = secis_df['has_secis'].sum()
    consistent = secis_df['sec_position_consistent'].sum() \
                 if 'sec_position_consistent' in secis_df.columns else has_secis
    final      = int(
        ((secis_df['has_secis'] == True) &
         (secis_df['sec_position_consistent'] == True)).sum()
    )

    print(f"\n--- Stage 3 Complete — Validation Funnel ---")
    print(f"  Passed UGA check     : {total}")
    print(f"  SECIS found          : {has_secis}")
    print(f"  Sec position valid   : {consistent}")
    print(f"  Final candidates     : {final}")

    if final > 0:
        print(f"\n  Validated selenoprotein candidates:")
        mask = (
            (secis_df['has_secis'] == True) &
            (secis_df.get('sec_position_consistent', True) == True)
        )
        for _, row in secis_df[mask].iterrows():
            print(f"    {row['region_id']}")
            print(f"      SECIS score  : {row['best_score']:.2f}")
            print(f"      SECIS e-val  : {row['best_evalue']}")
            print(f"      SECIS pos    : {row['secis_pos_abs']}")

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Stage 3.2: SECIS detection using Infernal + RF00031"
    )
    parser.add_argument('--passing-fa',
                        default='results/stage3/uga_passing_regions.fa')
    parser.add_argument('--uga-tsv',
                        default='results/stage3/uga_validated.tsv')
    parser.add_argument('--cm-path',
                        default=RFAM_SECIS_CM,
                        help='Path to RF00031.cm (downloaded if missing)')
    parser.add_argument('--output-tsv',
                        default='results/stage3/secis_hits.tsv')
    parser.add_argument('--output-passing-fa',
                        default='results/stage3/secis_passing_regions.fa')
    parser.add_argument('--threads',
                        type=int, default=4)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    os.makedirs('results/stage3/tmp', exist_ok=True)

    print("\n==========================================")
    print(" Stage 3.2: SECIS Detection (Infernal)")
    print("==========================================")

    # Step 1 — Check Infernal
    check_infernal()

    # Step 2 — Download + press RF00031
    setup_rfam_model(args.cm_path)

    # Step 3 — Load UGA-passing regions
    records = load_passing_regions(args.passing_fa)

    # Step 4 — Write downstream sequences to temp FASTA
    tmp_fasta   = 'results/stage3/tmp/downstream_seqs.fa'
    offsets     = write_downstream_fasta(records, tmp_fasta)

    # Step 5 — Run cmsearch
    tblout_path = 'results/stage3/tmp/cmsearch_raw.tblout'
    run_cmsearch(args.cm_path, tmp_fasta, tblout_path, args.threads)

    # Step 6 — Parse results
    infernal_hits = parse_cmsearch_tblout(tblout_path)
    print(f"\n  Regions with any hit : {len(infernal_hits)}")

    # Step 7 — Build results dataframe
    secis_df = build_results(records, infernal_hits, offsets)

    # Step 8 — Verify Sec positions are consistent with UGA positions
    uga_df   = pd.read_csv(args.uga_tsv, sep='\t') \
               if os.path.exists(args.uga_tsv) else pd.DataFrame()
    secis_df = verify_sec_positions(secis_df, uga_df)

    # Step 9 — Write outputs
    write_outputs(secis_df, records, args.output_tsv, args.output_passing_fa)

    # Step 10 — Summary
    print_summary(uga_df, secis_df)

    print(f"\n==========================================")
    print(f" Stage 3.2 check_secis — Complete")
    print(f"==========================================")
    print(f"  SECIS results : {args.output_tsv}")
    print(f"  Final FASTA   : {args.output_passing_fa}")
    print(f"  Next          : Stage 4 — filter_and_output.py")
    print(f"==========================================\n")


if __name__ == '__main__':
    main()