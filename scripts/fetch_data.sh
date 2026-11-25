#!/usr/bin/env bash
# Download the truth table, build the reference subset and T1K index, pick the samples and pull
# MHC reads for all of them. Six samples are sliced at a time; each takes a few minutes over HTTP.
# Usage: scripts/fetch_data.sh [JOBS]
set -euo pipefail
cd "$(dirname "$0")/.."
JOBS=${1:-6}
mkdir -p data/truth data/reads

TRUTH=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt
[[ -s data/truth/1000G_HLA_types.txt ]] || curl -fsSL "$TRUTH" -o data/truth/1000G_HLA_types.txt

scripts/fetch_reference.sh data
[[ -s assets/samples.tsv ]] || python scripts/select_samples.py

# HG001_GIAB is the 300x GIAB BAM: subsample to ~30x (seed 7, fraction 0.10)
awk -F'\t' 'NR > 1 {print $1, $6, ($1 == "HG001_GIAB" ? "7.10" : "")}' assets/samples.tsv \
  | xargs -P"$JOBS" -L1 sh -c './scripts/extract_hla_reads.sh "$0" "$1" data/reads data/ref/GRCh38_chr6_nonprimary.fa $2'
rm -f ./*.crai ./*.bai     # htslib caches remote indexes in the working directory

{ echo "sample,fastq_1,fastq_2"
  awk -F'\t' 'NR > 1 {print $1 "," $1 "_R1.fastq.gz," $1 "_R2.fastq.gz"}' assets/samples.tsv
} > data/reads/samplesheet.csv
echo "reads ready: $(ls data/reads/*_R2.fastq.gz | wc -l) samples"
