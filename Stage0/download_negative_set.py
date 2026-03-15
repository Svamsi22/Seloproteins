# download_negative_set.py
# Downloads Cys-homologs for each selenoprotein family

import requests
import argparse
import time
import os
from Bio import SeqIO
from io import StringIO

# Known selenoprotein families and their Cys-homolog search terms
# These are the well-characterized families where U/C exchange is documented
SELENO_FAMILIES = {
    "GPx":    "Glutathione peroxidase",
    "TXNRD":  "Thioredoxin reductase",
    "MSRB":   "Methionine sulfoxide reductase B",
    "DIO":    "Iodothyronine deiodinase",
    "SELENOW": "Selenoprotein W",
    "SELENOT": "Selenoprotein T",
    "SELENOH": "Selenoprotein H",
}

def download_cys_homologs_for_family(family_name, search_term):
    """
    For a given selenoprotein family, find the Cys-containing homologs.
    These are proteins from the same family but WITHOUT selenocysteine.
    """
    print(f"  Fetching Cys-homologs for {family_name}...")

    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "fasta",
        "query": (
            f"(protein_name:{search_term}) "
            f"AND (reviewed:true) "
            f"NOT (keyword:KW-0712)"   # explicitly exclude selenoproteins
        ),
        "size": "50"   # limit per family — we don't need thousands
    }

    try:
        response = requests.get(url, params=params, timeout=60)
        response.raise_for_status()
        records = list(SeqIO.parse(StringIO(response.text), "fasta"))

        # Tag each record with its family
        for r in records:
            r.description += f" [family={family_name}][negative_set]"

        print(f"    Found {len(records)} Cys-homologs")
        time.sleep(0.5)   # be polite to UniProt API
        return records

    except Exception as e:
        print(f"    WARNING: Failed for {family_name}: {e}")
        return []

def validate_negative_set(records):
    """
    Verify negative set does NOT contain U.
    If any slipped through (rare), discard them.
    """
    clean     = []
    has_U     = []
    MIN_LEN = 50

    for record in records:
        seq = str(record.seq).upper()
        if 'U' in seq:
            has_U.append(record)    # this shouldn't happen but check anyway 
        elif len(seq) >= MIN_LEN:
            clean.append(record)

    print(f"\n  Negative set validation:")
    print(f"    Clean (no U)     : {len(clean)}  ← kept")
    print(f"    Has U (rejected) : {len(has_U)}  ← discarded (are selenoproteins!)")

    return clean

def print_negative_summary(records):
    """Show what families are represented"""
    from collections import Counter

    families = []
    for r in records:
        desc = r.description
        if '[family=' in desc:
            fam = desc.split('[family=')[1].split(']')[0]
            families.append(fam)

    counts = Counter(families)
    print("\n  Negative set composition:")
    for fam, count in counts.most_common():
        print(f"    {fam:<12} : {count} Cys-homologs")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='data/processed/negative_set.fa')
    args = parser.parse_args()
    os.makedirs("data/processed", exist_ok=True)
    

    print("\n==========================================")
    print(" Building Negative Set (Cys-homologs)")
    print("==========================================")

    all_records = []
    for family_name, search_term in SELENO_FAMILIES.items():
        records = download_cys_homologs_for_family(family_name, search_term)
        all_records.extend(records)

    clean_records = validate_negative_set(all_records)
    print_negative_summary(clean_records)

    SeqIO.write(clean_records, args.output, "fasta")

    print(f"\n  Negative set saved : {args.output}")
    print(f"  Total sequences    : {len(clean_records)}")

if __name__ == '__main__':
    main()