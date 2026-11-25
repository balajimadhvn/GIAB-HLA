#!/usr/bin/env bash
# Reference resources:
#   * chr6 plus every non-primary contig (alts, decoys, HLA-* alleles, unplaced) of
#     GRCh38_full_analysis_set_plus_decoy_hla, the reference the NYGC 1000G CRAMs were aligned to. CRAM
#     containers on the HLA contigs also hold reads on neighbouring alt/decoy contigs, and htslib needs all
#     of them to decode a container. Fetched with ranged HTTP requests through the remote .fai (~300 Mb),
#     so the 3.1 GB genome is never downloaded.
#   * T1K HLA index built from the current IPD-IMGT/HLA release
# Needs samtools >= 1.10 and t1k-build.pl on PATH.
set -euo pipefail
OUT=${1:-data}
mkdir -p "$OUT/ref"
REF=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
FA=$OUT/ref/GRCh38_chr6_nonprimary.fa

curl -fsSL "$REF.fai" | awk '$1 == "chr6" || $1 !~ /^chr([0-9]+|X|Y|M)$/ {print $1}' > "$OUT/ref/contigs.txt"
touch "$FA"
[[ -s $FA.fai ]] || { [[ -s $FA ]] && samtools faidx "$FA"; } || true
missing=$(cut -f1 "$FA.fai" 2>/dev/null | grep -vxFf - "$OUT/ref/contigs.txt" || true)
if [[ -n $missing ]]; then
  echo "$missing" > "$OUT/ref/missing.txt"
  samtools faidx -r "$OUT/ref/missing.txt" "$REF" >> "$FA"
  rm -f "$OUT/ref/missing.txt" ./*.fa.fai
  samtools faidx "$FA"
fi
echo "reference subset: $(wc -l < "$FA.fai") contigs, $(du -h "$FA" | cut -f1)"

if [[ ! -s $OUT/ref/t1k_hla/t1k_hla_dna_seq.fa ]]; then
  # run t1k-build.pl with the perl from its own conda environment (it needs IO::Uncompress::Unzip)
  builder=$(command -v t1k-build.pl)
  ( cd "$OUT/ref" && "$(dirname "$builder")/perl" "$builder" -o t1k_hla --download IPD-IMGT/HLA > t1k_build.log 2>&1 )
  rm -f "$OUT/ref/t1k_hla/hla.dat.zip"
  ls "$OUT/ref/t1k_hla"
fi
