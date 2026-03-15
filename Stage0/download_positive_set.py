# download_positive_set.py
# Downloads selenoproteins from UniProt + SelenoDB 2.0

import requests
import argparse
import time
import os
from Bio import SeqIO
from io import StringIO
def download_uniprot_selenoproteins():
    """
    Download reviewed selenoproteins from UniProt Swiss-Prot.
    Uses two filters combined:
      KW-0712 = Selenoprotein keyword
      ft_non_std = SEC = has actual Sec feature evidence
    """
    print("\n[UniProt] Downloading reviewed selenoproteins...")

    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "fasta",
        "query": (
            "keyword:KW-0712 AND reviewed:true"        # actual SEC evidence in sequence
        )
    }

    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()

    records = list(SeqIO.parse(StringIO(response.text), "fasta"))
    print(f"  UniProt: {len(records)} selenoproteins downloaded")
    return records

def download_selenodb():
    """
    Download from SelenoDB 2.0 — industry standard curated set.
    SelenoDB contains manually verified selenoprotein sequences
    used by selenoprofiles4 as reference profiles.
    """
    print("\n[SelenoDB 2.0] Downloading curated selenoproteins...")

    # SelenoDB 2.0 human selenoproteome FASTA
    url = "https://selenodb.org/download/selenodb2.0_human.fa"

    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        records = list(SeqIO.parse(StringIO(response.text), "fasta"))
        print(f"  SelenoDB: {len(records)} selenoproteins downloaded")
        return records
    except Exception as e:
        print(f"  WARNING: SelenoDB download failed: {e}")
        print("  Continuing with UniProt only...")
        return []

def merge_and_deduplicate(uniprot_records, selenodb_records):
    """
    Merge both sources and remove duplicates by sequence identity.
    Keep track of source for each record.
    """
    seen_ids  = set()
    seen_seqs = set()
    merged    = []

    for record in uniprot_records:
        record.description += " [source=UniProt]"
        seq = str(record.seq).upper()
        if record.id not in seen_ids and seq not in seen_seqs:
            seen_ids.add(record.id)
            seen_seqs.add(seq)
            merged.append(record)

    added_from_selenodb = 0
    for record in selenodb_records:
        record.description += " [source=SelenoDB2.0]"
        seq = str(record.seq).upper()
        if record.id not in seen_ids and seq not in seen_seqs:
            seen_ids.add(record.id)
            seen_seqs.add(seq)
            merged.append(record)
            added_from_selenodb += 1

    print(f"\n  Merged positive set:")
    print(f"    From UniProt    : {len(uniprot_records)}")
    print(f"    Added from SelenoDB (new): {added_from_selenodb}")
    print(f"    Total unique    : {len(merged)}")
    return merged

def validate_positive_set(records):
    """
    Verify every sequence in positive set actually has U.
    Flag any that don't — these are annotation errors.
    """
    valid   = []
    invalid = []
    MIN_LEN = 50

    for record in records:
        seq = str(record.seq).upper()
        if 'U' in seq and len(seq) >= MIN_LEN:
            u_positions = [i+1 for i, aa in enumerate(seq) if aa == 'U']
            valid.append((record, u_positions))
        else:
            invalid.append(record)

    print(f"\n  Positive set validation:")
    print(f"    Valid (have U)  : {len(valid)}")
    print(f"    Invalid (no U)  : {len(invalid)}  ← discarded")

    if invalid:
        print("  Discarded (no U found in sequence):")
        for r in invalid:
            print(f"    {r.id}")

    return [r for r, _ in valid]

def main():
    print("SCRIPT STARTED")
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='data/processed/positive_set.fa')
    args = parser.parse_args()
    os.makedirs("data/processed", exist_ok=True) 

    uniprot_records  = download_uniprot_selenoproteins()
    selenodb_records = download_selenodb()

    merged   = merge_and_deduplicate(uniprot_records, selenodb_records)
    filtered = validate_positive_set(merged)

    SeqIO.write(filtered, args.output, "fasta")
    print(f"\n  Positive set saved: {args.output}")
    print(f"  Total sequences   : {len(filtered)}")

if __name__ == '__main__':
    main()