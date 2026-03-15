
# run_mmseqs2.py
# Stage 1: MMseqs2 translated search — genome vs selenoprotein database


import argparse
import subprocess
import shutil
import sys
import os
import time
import pandas as pd
from pathlib import Path


MMSEQS_COLUMNS = [
    "query_id",
    "target_protein",
    "identity",
    "alignment_length",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "evalue",
    "bitscore",
    "query_length",
    "target_length",
    "offset"
]

MMSEQS_FORMAT_OUTPUT = (
    "query,target,pident,alnlen,"
    "qstart,qend,tstart,tend,"
    "evalue,bits,qlen,tlen"
)


# ─────────────────────────────────────────────
# CHECK TOOL
# ─────────────────────────────────────────────

def check_mmseqs2():

    if not shutil.which("mmseqs"):
        print("ERROR: MMseqs2 not installed")
        sys.exit(1)

    result = subprocess.run(
        ["mmseqs", "version"],
        capture_output=True,
        text=True
    )

    print("MMseqs2 version:", result.stdout.strip())


# ─────────────────────────────────────────────
# SPLIT GENOME
# ─────────────────────────────────────────────

def split_genome(genome, chunk_dir, chunk_size=5_000_000):

    print("\n--- Splitting genome ---")

    Path(chunk_dir).mkdir(parents=True, exist_ok=True)

    header = None
    seq = ""

    with open(genome) as f:

        for line in f:

            if line.startswith(">"):
                header = line.strip()[1:]
            else:
                seq += line.strip()

    chunks = []
    length = len(seq)

    for i in range(0, length, chunk_size):

        chunk_seq = seq[i:i+chunk_size]

        chunk_file = f"{chunk_dir}/chunk_{i}.fa"

        with open(chunk_file, "w") as out:

            out.write(f">{header}_{i}\n")

            for j in range(0, len(chunk_seq), 60):
                out.write(chunk_seq[j:j+60] + "\n")

        chunks.append(chunk_file)

    print("Genome length:", length)
    print("Chunks:", len(chunks))

    return chunks


# ─────────────────────────────────────────────
# RUN SEARCH ON ONE CHUNK
# ─────────────────────────────────────────────

def run_mmseqs_chunk(query_fasta, target_db, work_dir, threads, sensitivity):

    offset = int(Path(query_fasta).stem.split("_")[1])

    query_db = f"{work_dir}/queryDB"
    result_db = f"{work_dir}/resultDB"
    tmp_dir = f"{work_dir}/tmp"

    subprocess.run(["mmseqs", "createdb", query_fasta, query_db], check=True)

    sensitivity_map = {
        "default": "4.0",
        "more-sensitive": "4.0",
        "ultra-sensitive": "5.7"
    }

    s = sensitivity_map[sensitivity]

    try:
        
        subprocess.run([
            "mmseqs","search",
            query_db,
            target_db,
            result_db,
            tmp_dir,
            "--search-type","2",
            "-s",s,
            "--threads",str(threads),
            "-e","1e-10",
            "--max-seqs","5",
            "--min-seq-id","0.5"
        ], check=True)
    except subprocess.CalledProcessError:
        print("No ORFs found in chunk — skipping")
        return None

    raw_tsv = f"{work_dir}/raw.tsv"

    subprocess.run([
        "mmseqs", "convertalis",
        query_db,
        target_db,
        result_db,
        raw_tsv,
        "--format-output",
        MMSEQS_FORMAT_OUTPUT
    ], check=True)

    # add offset column
    if os.path.exists(raw_tsv) and os.path.getsize(raw_tsv) > 0:

        df = pd.read_csv(raw_tsv, sep="\t", header=None)
        df["offset"] = offset
        df.to_csv(raw_tsv, sep="\t", header=False, index=False)

    return raw_tsv


# ─────────────────────────────────────────────
# PARSE OUTPUT
# ─────────────────────────────────────────────

def parse_and_clean(output_raw, output_tsv):

    if not os.path.exists(output_raw) or os.path.getsize(output_raw) == 0:

        pd.DataFrame(
            columns=MMSEQS_COLUMNS + ["coverage", "strand", "source"]
        ).to_csv(output_tsv, sep="\t", index=False)

        return 0

    df = pd.read_csv(output_raw, sep="\t", header=None)
    df.columns = MMSEQS_COLUMNS

    # Fix genome coordinates
    df["query_start"] = df["query_start"] + df["offset"]
    df["query_end"]   = df["query_end"] + df["offset"]

    df["aligned_protein_len"] = (
        df["target_end"] - df["target_start"] + 1
    ).abs()

    df["coverage"] = (
        df["aligned_protein_len"] /
        df["target_length"] * 100
    ).round(2)

    df["strand"] = df.apply(
        lambda r: "-" if r["query_start"] > r["query_end"] else "+",
        axis=1
    )

    df["query_start"], df["query_end"] = (
        df[["query_start", "query_end"]].min(axis=1),
        df[["query_start", "query_end"]].max(axis=1)
    )

    if df["identity"].max() <= 1:
        df["identity"] = (df["identity"] * 100).round(2)

    df["source"] = "MMseqs2"

    df.drop(columns=["offset"], inplace=True)

    df.to_csv(output_tsv, sep="\t", index=False)

    return len(df)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genome",
        default="data/raw/chr22_clean.fa"
    )

    parser.add_argument(
        "--target-db",
        default="databases/mmseqs_db/selenoDB"
    )

    parser.add_argument(
        "--output-raw",
        default="results/stage1/mmseqs_raw.tsv"
    )

    parser.add_argument(
        "--output-tsv",
        default="results/stage1/hits_mmseqs2.tsv"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=2
    )

    parser.add_argument(
        "--sensitivity",
        default="more-sensitive",
        choices=["default", "more-sensitive", "ultra-sensitive"]
    )

    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_raw), exist_ok=True)

    print("\n==========================================")
    print(" Stage 1: MMseqs2 Search")
    print("==========================================")

    check_mmseqs2()

    chunks = split_genome(args.genome, "tmp_chunks")

    start = time.time()

    outputs = []

    for i, chunk in enumerate(chunks):

        print(f"\nRunning MMseqs2 on chunk {i+1}/{len(chunks)}")

        work = f"tmp_mmseqs_{i}"

        os.makedirs(work, exist_ok=True)

        raw = run_mmseqs_chunk(
            chunk,
            args.target_db,
            work,
            args.threads,
            args.sensitivity
        )

        if raw:
            outputs.append(raw)

    print("\n--- Merging results ---")

    with open(args.output_raw, "w") as out:

        for f in outputs:

            if os.path.exists(f):

                with open(f) as r:
                    out.write(r.read())

    elapsed = time.time() - start

    hits = parse_and_clean(args.output_raw, args.output_tsv)

    print("\n==========================================")
    print(" Stage 1 MMseqs2 — Complete")
    print("==========================================")
    print("Hits:", hits)
    print("Runtime:", round(elapsed,1),"s")
    print("Output:", args.output_tsv)
    print("==========================================\n")


if __name__ == "__main__":
    main()