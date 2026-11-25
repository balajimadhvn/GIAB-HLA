#!/usr/bin/env python3
"""
Choose the benchmark samples and write assets/samples.tsv.

  * NA12878 (GIAB HG001) and both CEPH parents, NA12891 and NA12892, for a trio consistency check
  * HG001 again from the GIAB 300x HiSeq data, downsampled to 30x, to compare two independent data sets
  * 4 more samples per 1000 Genomes super-population (AFR, AMR, EAS, EUR, SAS), drawn with a fixed
    hash so the choice is reproducible, restricted to samples whose 10 truth alleles are all typed

Truth: Abi-Rached et al. 2018 (PLoS One 13:e0206512), HLA types for the 1000 Genomes panel.
Reads: NYGC 30x GRCh38 CRAMs (1000 Genomes high-coverage release, AWS Open Data mirror).
"""
import csv
import hashlib
import urllib.request
from pathlib import Path

BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
INDEXES = ["1000G_2504_high_coverage.sequence.index", "1000G_698_related_high_coverage.sequence.index"]
S3 = "https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/"
GIAB_HG001 = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/"
              "NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam")
TRIO = ["NA12878", "NA12891", "NA12892"]
PER_REGION = 4


def cram_urls() -> dict[str, str]:
    urls = {}
    for idx in INDEXES:
        related = "698" in idx
        with urllib.request.urlopen(BASE + idx) as fh:
            for line in fh.read().decode().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                ftp = line.split("\t")[0]
                run, name = ftp.split("/")[-2], ftp.split("/")[-1]
                sample = name.split(".")[0]
                sub = "additional_698_related/data/" if related else "data/"
                urls[sample] = f"{S3}{sub}{run}/{name}"
    return urls


def main() -> None:
    truth = list(csv.DictReader(open("data/truth/1000G_HLA_types.txt"), delimiter="\t"))
    allele_cols = [c for c in truth[0] if c.startswith("HLA-")]
    complete = [r for r in truth if all(r[c].strip() not in ("", "None") for c in allele_cols)]
    urls = cram_urls()

    rows = []
    for s in TRIO:
        r = next(x for x in truth if x["Sample ID"] == s)
        rows.append(dict(sample=s, region=r["Region"], population=r["Population"], truth_id=s,
                         source="NYGC 30x", url=urls[s], group="trio"))
    rows.append(dict(sample="HG001_GIAB", region="EUR", population="CEU", truth_id="NA12878",
                     source="GIAB 300x HiSeq, downsampled to 30x", url=GIAB_HG001, group="trio"))

    key = lambda r: hashlib.sha1(("giab-hla-v1:" + r["Sample ID"]).encode()).hexdigest()
    for region in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        pool = sorted((r for r in complete if r["Region"] == region and r["Sample ID"] in urls
                       and r["Sample ID"] not in TRIO), key=key)
        for r in pool[:PER_REGION]:
            rows.append(dict(sample=r["Sample ID"], region=region, population=r["Population"],
                             truth_id=r["Sample ID"], source="NYGC 30x", url=urls[r["Sample ID"]],
                             group="population"))

    out = Path("assets/samples.tsv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}: {len(rows)} samples")


if __name__ == "__main__":
    main()
