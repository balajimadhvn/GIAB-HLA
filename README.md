# GIAB-HLA

[![CI](https://github.com/balajimadhvn/GIAB-HLA/actions/workflows/ci.yml/badge.svg)](https://github.com/balajimadhvn/GIAB-HLA/actions/workflows/ci.yml)

HLA class I and II typing from 30x whole-genome short reads, using T1K and OptiType in a Nextflow pipeline.
It is benchmarked against published HLA types for GIAB HG001 (NA12878), the CEPH trio, and 20 more
1000 Genomes samples from all five super-populations.

It is the HLA companion to [GermlineQC](https://github.com/balajimadhvn/GermlineQC). That repo covers the
germline small-variant pipeline, its GIAB validation and ACMG classification.

**Report:** https://balajimadhvn.github.io/GIAB-HLA/

```
remote CRAM/BAM ── samtools -M -L (MHC + HLA/alt contigs) ── paired FASTQ
                                                               ├─ T1K      HLA-A, -B, -C, -DQB1, -DRB1
                                                               └─ OptiType HLA-A, -B, -C
                                                                     │
                         2-field calls ── concordance vs truth, tool agreement, trio check ── R report
```

## Results

24 genomes: 23 people, with HG001 typed twice from independent data sets. Each locus pairs the two calls with the
two truth alleles in whichever order gives more matches, and the score is the fraction of alleles that match.
Intervals are 95% Wilson intervals.

| Tool | Loci | 1-field | 2-field |
|---|---|---|---|
| T1K | A, B, C, DQB1, DRB1 | 239/240 (99.6%, 97.7-99.9) | 237/240 (98.8%, 96.4-99.6) |
| OptiType | A, B, C | 144/144 (100%, 97.4-100) | 143/144 (99.3%, 96.2-99.9) |

| Locus | T1K 2-field | OptiType 2-field |
|---|---|---|
| HLA-A | 47/48 | 47/48 |
| HLA-B | 48/48 | 48/48 |
| HLA-C | 48/48 | 48/48 |
| HLA-DQB1 | 47/48 | - |
| HLA-DRB1 | 47/48 | - |

- **Tool agreement:** T1K and OptiType give the same 2-field class I genotype for all 72 sample-locus pairs.
- **Trio:** NA12878's genotype can be built from one allele of NA12891 and one of NA12892 at every locus, for both
  tools. This check does not depend on the truth set.
- **The three misses:**
  - NA18618 HLA-A: both tools call A\*24:02 where the truth has A\*24:01. Two independent methods agree, so the
    truth entry is the more likely error. It was typed from exome data against an older IMGT/HLA release.
  - HG001_GIAB HLA-DQB1: T1K calls 05:01/05:03 against 02:01/05:01. The NYGC 30x data for the same person
    (NA12878) is typed correctly, so this is a data-set effect: different library, aligner and depth after
    downsampling.
  - HG01631 HLA-DRB1: T1K calls 04:14 against 04:02. The allele group is right; the second field is wrong.

Tables: [`docs/tables/`](docs/tables). They hold accuracy by locus, per-sample concordance, tool agreement, the
trio check and every call.

### Runtime

All runs were on a 16 GB laptop, with the Nextflow executor capped at 8 CPUs and 6.5 GB.

| Step | Per sample | Peak memory |
|---|---|---|
| Read extraction (remote CRAM, about 1M reads) | 1-3 min (HG001 300x BAM: about 15 min) | small |
| T1K | about 4 min | 3.2 GB |
| OptiType | 2-4 min | 1.3 GB |

## How the reads are obtained

The NYGC 1000 Genomes CRAMs are not downloaded. `scripts/extract_hla_reads.sh` makes HTTP range requests for
chr6:29.6-33.5 Mb plus every HLA-\* and chr6 alt contig. It uses samtools' multi-region iterator (`-M -L`), so
the ~540 regions are read in one pass over the index rather than one seek each. Decoding a CRAM slice still
needs the reference for every contig in a container, and containers mix contigs. For that reason
`scripts/fetch_reference.sh` pulls chr6 and all non-primary contigs (3,342 sequences, about 290 MB) from the
remote FASTA index, not the full 3.1 GB genome.

HG001 also comes from the GIAB 300x HiSeq BAM, downsampled to about 30x (`samtools view -s 7.10`). That gives
a second, independent data set for one person.

## Usage

```bash
mamba env create -f environment.yml            # Nextflow, samtools, Python, R
mamba env create -f envs/t1k.yml
mamba env create -f envs/optitype.yml
conda activate giab-hla

scripts/fetch_data.sh                          # truth table, reference subset, T1K index, reads
nextflow run main.nf                           # results/ and results/hla_report.html
```

To use the existing environments rather than letting Nextflow build them, pass
`--t1k_env`, `--optitype_env` and `--base_env` with the environment paths.

Smoke test: 5,831 real HG001 read pairs and a reduced T1K index, both shipped in `tests/data`. The reduced index holds one allele per 2-field group, so the test checks that the pipeline
runs, not typing accuracy. OptiType is still expected to get 6/6 on it, and CI checks that.

```bash
nextflow run main.nf -profile test
pytest
```

### Your own samples

`--input` takes a CSV with `sample,fastq_1,fastq_2` columns. Paths are relative to the CSV. `--samples` and
`--truth` are only needed for scoring. Without a truth row a sample is still typed, and its calls are in
`results/tables/calls.tsv`.

## Repository

| Path | Contents |
|---|---|
| `main.nf`, `nextflow.config` | T1K, OptiType, summary and report processes; `test` profile |
| `bin/hla_concordance.py` | Allele normalisation, truth ambiguity handling, accuracy, agreement, trio check |
| `reports/hla_report.Rmd` | R Markdown report (ggplot2) |
| `scripts/` | Sample selection, reference subset, remote read extraction, data download |
| `assets/samples.tsv` | The 24 benchmark samples, populations and source URLs |
| `envs/`, `environment.yml` | Pinned conda environments |
| `tests/` | pytest unit tests and the smoke-test data |
| `docs/` | Published report and tables |

## Limitations

- The truth set was itself made computationally, from exome data and IMGT/HLA 3.28. It is a reference, not
  a gold standard; NA18618 HLA-A is probably an example.
- 48 alleles per locus give 95% intervals about 10 points wide. A clinical HLA validation would need more
  samples, sequence-based typing as the truth, and coverage of rare alleles.
- HLA-DPB1, -DQA1 and -DRB3/4/5 are not scored because the truth set does not include them. T1K still
  reports them.
- Short-read typing is limited at 3rd and 4th field resolution and for novel alleles.

## Data sources

- HLA truth: Abi-Rached et al. 2018, *PLoS One* 13:e0206512 (1000 Genomes HLA types, 20181129 release)
- Reads: NYGC 30x GRCh38 1000 Genomes high-coverage CRAMs (AWS Open Data); GIAB HG001 300x HiSeq (NCBI)
- Reference: GRCh38_full_analysis_set_plus_decoy_hla (1000 Genomes)
- HLA alleles: IPD-IMGT/HLA, current release at build time
- Tools: [T1K](https://github.com/mourisl/T1K) 1.0.10, [OptiType](https://github.com/FRED-2/OptiType) 1.3.5

## License

MIT
