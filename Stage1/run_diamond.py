
# run_diamond.py
# Stage 1: DIAMOND blastx search — genome vs selenoprotein database
# Fixed version with correct coordinate offsets

import argparse
import subprocess
import shutil
import sys
import os
import time
import pandas as pd
from pathlib import Path

# ─────────────────────────────────────────────
# DIAMOND OUTPUT FORMAT
# ─────────────────────────────────────────────

DIAMOND_COLUMNS = [
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
    "offset",
]

DIAMOND_OUTFMT = (
    "6 "
    "qseqid sseqid pident length "
    "qstart qend sstart send "
    "evalue bitscore qlen slen"
)

# ─────────────────────────────────────────────
# CHECK DIAMOND
# ─────────────────────────────────────────────

def check_diamond():
    if not shutil.which("diamond"):
        print("ERROR: DIAMOND not installed.")
        sys.exit(1)

    result = subprocess.run(
        ["diamond", "version"],
        capture_output=True,
        text=True
    )
    print("  DIAMOND version:", result.stdout.split("\n")[0])


# ─────────────────────────────────────────────
# SPLIT GENOME INTO CHUNKS
# ─────────────────────────────────────────────

def split_genome(genome, chunk_dir, chunk_size=5_000_000):

    print("\n--- Splitting genome into chunks ---")

    Path(chunk_dir).mkdir(parents=True, exist_ok=True)

    seq = ""
    header = None
    chunks = []

    with open(genome) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
            else:
                seq += line.strip()

    length = len(seq)

    for i in range(0, length, chunk_size):
        chunk_seq = seq[i:i + chunk_size]
        chunk_file = f"{chunk_dir}/chunk_{i}.fa"

        with open(chunk_file, "w") as out:
            out.write(f">{header[1:]}_{i}\n")
            for j in range(0, len(chunk_seq), 60):
                out.write(chunk_seq[j:j+60] + "\n")

        chunks.append(chunk_file)

    print(f"  Genome length : {length:,}")
    print(f"  Chunks created: {len(chunks)}")

    return chunks


# ─────────────────────────────────────────────
# RUN DIAMOND ON ONE CHUNK
# ─────────────────────────────────────────────

def run_diamond_chunk(chunk, db, output, threads, sensitivity):

    offset = int(Path(chunk).stem.split("_")[1])

    cmd = [
        "diamond", "blastx",
        "--db", db,
        "--query", chunk,
        "--out", output,
        "--outfmt", *DIAMOND_OUTFMT.split(),
        "--threads", str(threads),
        f"--{sensitivity}",
        "--evalue", "1e-5"
    ]

    subprocess.run(cmd, check=True)

    # add offset column to output
    if os.path.exists(output) and os.path.getsize(output) > 0:

        df = pd.read_csv(output, sep="\t", header=None)
        df["offset"] = offset
        df.to_csv(output, sep="\t", header=False, index=False)


# ─────────────────────────────────────────────
# RUN DIAMOND ON ALL CHUNKS
# ─────────────────────────────────────────────

def run_diamond(genome, db, output_raw, threads, sensitivity):

    chunk_dir = "tmp_chunks"
    out_dir = "tmp_hits"

    Path(out_dir).mkdir(exist_ok=True)

    chunks = split_genome(genome, chunk_dir)

    start = time.time()

    outputs = []

    for i, chunk in enumerate(chunks):

        out_file = f"{out_dir}/hits_{i}.tsv"

        print(f"\nRunning DIAMOND on chunk {i+1}/{len(chunks)}")

        run_diamond_chunk(
            chunk,
            db,
            out_file,
            threads,
            sensitivity
        )

        outputs.append(out_file)

    print("\n--- Merging results ---")

    with open(output_raw, "w") as outfile:
        for f in outputs:
            if os.path.exists(f) and os.path.getsize(f) > 0:
                with open(f) as infile:
                    outfile.write(infile.read())

    elapsed = time.time() - start

    return elapsed


# ─────────────────────────────────────────────
# PARSE OUTPUT AND FIX COORDINATES
# ─────────────────────────────────────────────

def parse_and_clean(output_raw, output_tsv):

    if not os.path.exists(output_raw) or os.path.getsize(output_raw) == 0:

        print("WARNING: no DIAMOND hits found")

        pd.DataFrame(
            columns=DIAMOND_COLUMNS + ["coverage", "strand", "source"]
        ).to_csv(output_tsv, sep="\t", index=False)

        return 0

    df = pd.read_csv(output_raw, sep="\t", header=None)
    df.columns = DIAMOND_COLUMNS

    # ── extract chunk offset from query id ──
    # example: chr22_10000000
    # df["offset"] = df["query_id"].astype(str).str.extract(r"_(\d+)$").astype(int)

    # ── fix coordinates ──
    df["query_start"] = df["query_start"] + df["offset"]
    df["query_end"] = df["query_end"] + df["offset"]

    # coverage calculation
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
        df[["query_start", "query_end"]].max(axis=1),
    )

    df["source"] = "DIAMOND"

    df.drop(columns=["offset"], inplace=True)

    df.to_csv(output_tsv, sep="\t", index=False)

    return len(df)


# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────

def print_summary(output_tsv, elapsed):

    if not os.path.exists(output_tsv):
        return

    df = pd.read_csv(output_tsv, sep="\t")

    if df.empty:
        print("\nNo hits found.")
        return

    print("\n--- DIAMOND Results Summary ---")
    print("  Total hits:", len(df))
    print("  Unique proteins:", df["target_protein"].nunique())
    print("  Mean identity:", round(df["identity"].mean(), 2))
    print("  Mean coverage:", round(df["coverage"].mean(), 2))
    print("  Best evalue:", df["evalue"].min())
    print("  Runtime:", round(elapsed, 1), "seconds")


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
        "--diamond-db",
        default="databases/diamond_db/seleno"
    )

    parser.add_argument(
        "--output-raw",
        default="results/stage1/diamond_raw.tsv"
    )

    parser.add_argument(
        "--output-tsv",
        default="results/stage1/hits_diamond.tsv"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=2
    )

    parser.add_argument(
        "--sensitivity",
        default="more-sensitive",
        choices=["sensitive", "more-sensitive", "ultra-sensitive"]
    )

    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output_raw), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    print("\n==========================================")
    print(" Stage 1: DIAMOND blastx Search")
    print("==========================================")

    check_diamond()

    elapsed = run_diamond(
        args.genome,
        args.diamond_db,
        args.output_raw,
        args.threads,
        args.sensitivity
    )

    n_hits = parse_and_clean(
        args.output_raw,
        args.output_tsv
    )

    print_summary(args.output_tsv, elapsed)

    print("\n==========================================")
    print(" Stage 1 DIAMOND — Complete")
    print("==========================================")
    print("Hits found:", n_hits)
    print("Output:", args.output_tsv)
    print("==========================================\n")


if __name__ == "__main__":
    main()