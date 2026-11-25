#!/usr/bin/env python3
"""
Compare HLA calls from T1K and OptiType with published HLA types.

Outputs (tab-separated, in --outdir):
  calls.tsv         every call, normalised to 2-field resolution (e.g. A*01:01)
  concordance.tsv   per sample, tool and locus: number of alleles matching the truth (0, 1 or 2)
  accuracy.tsv      allele-level accuracy per tool and locus, at 1-field and 2-field resolution
  agreement.tsv     T1K vs OptiType agreement for class I
  trio.tsv          Mendelian consistency of the NA12878 / NA12891 / NA12892 calls

Truth cells can hold ambiguity ("03:01/04" = 03:01 or 03:04) and a trailing '*'. Empty and "None" cells
are untyped and left out of the denominator.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from itertools import permutations
from pathlib import Path

LOCI = ["A", "B", "C", "DQB1", "DRB1"]
CLASS_I = ["A", "B", "C"]


# ------------------------------------------------------------------ allele handling

def two_field(allele: str) -> str | None:
    """'HLA-A*01:01:01:01' / 'A*01:01' / '01:01' -> '01:01'. Returns None for empty calls."""
    if not allele or allele in (".", "None", "NA"):
        return None
    a = allele.strip().split(",")[0].removeprefix("HLA-")   # T1K lists equally likely alleles, comma-separated
    a = a.split("*", 1)[1] if "*" in a else a
    a = re.sub(r"[A-Za-z]$", "", a)          # expression suffixes (N, L, Q ...)
    fields = a.split(":")
    return ":".join(fields[:2]) if len(fields) >= 2 else None


def truth_options(cell: str) -> set[str] | None:
    """'03:02/03/04' -> {'03:02', '03:03', '03:04'}; '02:02*' -> {'02:02'}."""
    cell = cell.strip().rstrip("*")
    if cell in ("", "None"):
        return None
    first, *rest = cell.split("/")
    group = first.split(":")[0]
    out = {first}
    for r in rest:
        out.add(r if ":" in r else f"{group}:{r}")
    return out


def match(call: str, options: set[str], level: int) -> bool:
    if level == 1:
        return call.split(":")[0] in {o.split(":")[0] for o in options}
    return call in options


def best_pairing(calls: list[str], truth: list[set[str]], level: int) -> int:
    """Number of correct alleles under the best assignment of 2 calls to 2 truth alleles."""
    if len(calls) != 2 or len(truth) != 2:
        return 0
    return max(sum(match(c, t, level) for c, t in zip(p, truth)) for p in permutations(calls))


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


# ------------------------------------------------------------------ tool parsers

def parse_t1k(path: Path) -> dict[str, list[str]]:
    """T1K *_genotype.tsv: gene, n_alleles, allele1, abundance1, quality1, allele2, abundance2, quality2."""
    calls = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        locus = f[0].removeprefix("HLA-")
        if locus not in LOCI:
            continue
        a1, a2 = two_field(f[2]), two_field(f[5])
        if a1 is None:
            continue
        calls[locus] = [a1, a2 if a2 else a1]       # one allele reported -> homozygous
    return calls


def parse_optitype(path: Path) -> dict[str, list[str]]:
    rows = list(csv.reader(path.open(), delimiter="\t"))
    if len(rows) < 2:
        return {}
    rec = dict(zip(rows[0], rows[1]))
    calls = {}
    for locus in CLASS_I:
        a1, a2 = two_field(rec.get(f"{locus}1", "")), two_field(rec.get(f"{locus}2", ""))
        if a1:
            calls[locus] = [a1, a2 or a1]
    return calls


def sample_from(path: Path, suffix: str) -> str:
    return path.name.removesuffix(suffix)


# ------------------------------------------------------------------ analyses

def read_truth(path: Path) -> dict[str, dict[str, list[set[str]] | None]]:
    truth = {}
    for r in csv.DictReader(path.open(), delimiter="\t"):
        loci = {}
        for locus in LOCI:
            a, b = truth_options(r.get(f"HLA-{locus} 1", "")), truth_options(r.get(f"HLA-{locus} 2", ""))
            loci[locus] = [a, b] if a and b else None
        truth[r["Sample ID"]] = loci
    return truth


def trio_consistent(child: list[str], father: list[str], mother: list[str]) -> bool:
    """Child carries one allele from each parent (2-field)."""
    return any(c1 in father and c2 in mother for c1, c2 in permutations(child))


def write(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1k", type=Path, nargs="+", required=True)
    ap.add_argument("--optitype", type=Path, nargs="+", required=True)
    ap.add_argument("--samples", type=Path, required=True)
    ap.add_argument("--truth", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    meta = {r["sample"]: r for r in csv.DictReader(a.samples.open(), delimiter="\t")}
    truth = read_truth(a.truth)
    tools = {
        "T1K": {sample_from(p, "_genotype.tsv"): parse_t1k(p) for p in a.t1k},
        "OptiType": {sample_from(p, "_result.tsv"): parse_optitype(p) for p in a.optitype},
    }

    calls_rows, conc_rows = [], []
    for tool, by_sample in tools.items():
        for sample, loci in sorted(by_sample.items()):
            m = meta.get(sample, {})
            t = truth.get(m.get("truth_id", sample), {})
            for locus, alleles in loci.items():
                calls_rows.append(dict(sample=sample, tool=tool, locus=locus,
                                       allele1=f"{locus}*{alleles[0]}", allele2=f"{locus}*{alleles[1]}"))
                tl = t.get(locus)
                if tl is None:
                    continue
                conc_rows.append(dict(sample=sample, region=m.get("region", ""), group=m.get("group", ""),
                                      tool=tool, locus=locus,
                                      called="/".join(alleles), truth="/".join("|".join(sorted(x)) for x in tl),
                                      match_1field=best_pairing(alleles, tl, 1),
                                      match_2field=best_pairing(alleles, tl, 2)))
            # loci with a truth but no call count as 0 correct
            for locus, tl in t.items():
                if tl is not None and locus not in loci and (tool == "T1K" or locus in CLASS_I):
                    conc_rows.append(dict(sample=sample, region=m.get("region", ""), group=m.get("group", ""),
                                          tool=tool, locus=locus, called="no call",
                                          truth="/".join("|".join(sorted(x)) for x in tl),
                                          match_1field=0, match_2field=0))
    write(a.outdir / "calls.tsv", calls_rows)
    write(a.outdir / "concordance.tsv", conc_rows)

    acc_rows = []
    for tool in tools:
        loci = LOCI if tool == "T1K" else CLASS_I
        for locus in loci + ["ALL"]:
            rows = [r for r in conc_rows if r["tool"] == tool and (locus == "ALL" or r["locus"] == locus)]
            n = 2 * len(rows)
            for level in (1, 2):
                k = sum(r[f"match_{level}field"] for r in rows)
                lo, hi = wilson(k, n)
                acc_rows.append(dict(tool=tool, locus=locus, resolution=f"{level}-field", correct=k, alleles=n,
                                     accuracy=k / n if n else math.nan, ci_lo=lo, ci_hi=hi))
    write(a.outdir / "accuracy.tsv", acc_rows)

    agree_rows = []
    for sample in sorted(set(tools["T1K"]) & set(tools["OptiType"])):
        for locus in CLASS_I:
            t1 = tools["T1K"][sample].get(locus)
            ot = tools["OptiType"][sample].get(locus)
            if t1 and ot:
                n_same = max(sum(x == y for x, y in zip(p, ot)) for p in permutations(t1))
                agree_rows.append(dict(sample=sample, locus=locus, t1k="/".join(t1), optitype="/".join(ot),
                                       shared_alleles=n_same))
    write(a.outdir / "agreement.tsv", agree_rows)

    trio_rows = []
    for tool, by_sample in tools.items():
        if not {"NA12878", "NA12891", "NA12892"} <= set(by_sample):
            continue
        child, father, mother = (by_sample[s] for s in ("NA12878", "NA12891", "NA12892"))
        for locus in sorted(set(child) & set(father) & set(mother), key=LOCI.index):
            trio_rows.append(dict(tool=tool, locus=locus, child="/".join(child[locus]),
                                  father="/".join(father[locus]), mother="/".join(mother[locus]),
                                  consistent=trio_consistent(child[locus], father[locus], mother[locus])))
    write(a.outdir / "trio.tsv", trio_rows)

    for r in acc_rows:
        if r["locus"] == "ALL":
            print(f"{r['tool']:<9} {r['resolution']}: {r['correct']}/{r['alleles']} = {r['accuracy']:.4f}")


if __name__ == "__main__":
    main()
