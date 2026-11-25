#!/usr/bin/env bash
# Pull reads for the MHC class I/II genes from a remote CRAM/BAM and write paired FASTQs.
# Regions: chr6 MHC (29.6-33.5 Mb, GRCh38) plus every HLA-* and chr6 alt contig in the reference,
# because the NYGC alignments place many HLA reads on the HLA allele contigs.
# Usage: extract_hla_reads.sh SAMPLE URL OUTDIR REF [SUBSAMPLE_FRACTION]
set -euo pipefail
sample=$1; url=$2; out=$3; ref=$4; frac=${5:-}
mkdir -p "$out"
[[ -s $out/${sample}_R2.fastq.gz ]] && exit 0

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
printf 'chr6\t29600000\t33500000\n' > "$tmp/regions.bed"
refopt=()
if [[ $url == *.cram ]]; then
  awk -v OFS='\t' '$1 ~ /^HLA-/ || $1 ~ /^chr6_.*_alt$/ {print $1, 0, $2}' "$ref.fai" >> "$tmp/regions.bed"
  refopt=(-T "$ref")
fi
sub=(); [[ -n $frac ]] && sub=(-s "$frac")
# -M: one multi-region pass over the index instead of a separate remote seek for each of the ~540 contigs
samtools view -M -L "$tmp/regions.bed" -@2 -b -F 0x900 "${sub[@]}" "${refopt[@]}" -o "$tmp/slice.bam" "$url"
samtools collate -@2 -u -O "$tmp/slice.bam" "$tmp/c" \
  | samtools fastq -@2 -s /dev/null -0 /dev/null -1 "$out/${sample}_R1.fastq.gz" -2 "$out/${sample}_R2.fastq.gz" - 2> "$out/${sample}.fastq.log"
echo "$sample $(grep -o 'processed [0-9]* reads' "$out/${sample}.fastq.log")"
