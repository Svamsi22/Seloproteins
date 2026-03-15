#!/bin/bash
# bin/download_genome.sh

set -euo pipefail

mkdir -p data/raw data/ground_truth

echo "=========================================="
echo " Stage 0: Genome + Ground Truth Download"
echo "=========================================="

# ── 1. Download chr22 FASTA ───────────────────────────────────
echo "[1/4] Downloading Human chr22 (GRCh38 release 111)..."
wget -q --show-progress -P data/raw/ \
  "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"

# ── 2. Decompress FASTA ───────────────────────────────────────
echo "[2/4] Decompressing chr22.fa..."
gunzip -f data/raw/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
mv data/raw/Homo_sapiens.GRCh38.dna.chromosome.22.fa data/raw/chr22.fa
echo "      Saved: data/raw/chr22.fa"

# ── 3. Download ground truth GFF3 ────────────────────────────
echo "[3/4] Downloading Ensembl chr22 GFF3 (ground truth)..."
wget -q --show-progress -P data/ground_truth/ \
  "https://ftp.ensembl.org/pub/release-111/gff3/homo_sapiens/Homo_sapiens.GRCh38.111.chromosome.22.gff3.gz"

# ── 4. Decompress GFF3 ───────────────────────────────────────
echo "[4/4] Decompressing GFF3..."
gunzip -f data/ground_truth/Homo_sapiens.GRCh38.111.chromosome.22.gff3.gz
mv data/ground_truth/Homo_sapiens.GRCh38.111.chromosome.22.gff3 \
   data/ground_truth/ensembl_chr22_full.gff3
echo "      Saved: data/ground_truth/ensembl_chr22_full.gff3"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "=========================================="
echo " Downloads Complete"
echo "=========================================="
echo "  data/raw/chr22.fa"
echo "  data/ground_truth/ensembl_chr22_full.gff3"
echo "=========================================="