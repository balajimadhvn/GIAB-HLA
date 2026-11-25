from pathlib import Path

import hla_concordance as h


def test_two_field_normalisation():
    assert h.two_field("HLA-A*01:01:01:01") == "01:01"
    assert h.two_field("A*24:02") == "24:02"
    assert h.two_field("HLA-DRB1*01:01:01,HLA-DRB1*01:01:38,HLA-DRB1*01:135") == "01:01"
    assert h.two_field("01:01") == "01:01"
    assert h.two_field("HLA-C*07:01:01:01N") == "07:01"
    assert h.two_field(".") is None and h.two_field("") is None


def test_truth_ambiguity_codes():
    assert h.truth_options("03:01/04") == {"03:01", "03:04"}
    assert h.truth_options("03:02/03/04") == {"03:02", "03:03", "03:04"}
    assert h.truth_options("02:02*") == {"02:02"}
    assert h.truth_options("None") is None and h.truth_options(" ") is None


def test_best_pairing_is_order_independent():
    truth = [{"01:01"}, {"11:01"}]
    assert h.best_pairing(["11:01", "01:01"], truth, 2) == 2
    assert h.best_pairing(["01:01", "01:01"], truth, 2) == 1
    assert h.best_pairing(["01:02", "11:04"], truth, 2) == 0
    assert h.best_pairing(["01:02", "11:04"], truth, 1) == 2


def test_homozygous_truth_needs_both_calls():
    truth = [{"05:01"}, {"05:01"}]
    assert h.best_pairing(["05:01", "05:01"], truth, 2) == 2
    assert h.best_pairing(["05:01", "06:02"], truth, 2) == 1


def test_trio_consistency():
    # NA12878 DRB1 01:01/03:01, father 03:01/15:01, mother 01:01/01:01
    assert h.trio_consistent(["01:01", "03:01"], ["03:01", "15:01"], ["01:01", "01:01"])
    assert not h.trio_consistent(["04:01", "03:01"], ["03:01", "15:01"], ["01:01", "01:01"])


def test_parsers(tmp_path):
    t1k = tmp_path / "S_genotype.tsv"
    t1k.write_text("HLA-A\t2\tHLA-A*01:01:01\t100\t60\tHLA-A*11:01:01\t90\t60\n"
                   "HLA-DRB1\t1\tHLA-DRB1*15:01:01\t80\t40\t.\t0\t-1\n"
                   "HLA-E\t1\tHLA-E*01:01:01\t50\t30\t.\t0\t-1\n")
    assert h.parse_t1k(t1k) == {"A": ["01:01", "11:01"], "DRB1": ["15:01", "15:01"]}
    opt = tmp_path / "S_result.tsv"
    opt.write_text("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
                   "0\tA*01:01\tA*11:01\tB*08:01\tB*56:01\tC*01:02\tC*07:01\t1000\t950.0\n")
    assert h.parse_optitype(opt)["B"] == ["08:01", "56:01"]


def test_wilson_bounds():
    lo, hi = h.wilson(48, 48)
    assert 0.9 < lo < 1.0 and abs(hi - 1.0) < 1e-9


def test_real_t1k_output():
    """T1K 1.0.10 output for HG001 (GIAB 300x BAM, MHC reads downsampled to ~30x)."""
    calls = h.parse_t1k(Path(__file__).parent / "data" / "HG001_GIAB_genotype.tsv")
    assert sorted(calls["A"]) == ["01:01", "11:01"]
    assert sorted(calls["DRB1"]) == ["01:01", "03:01"]
    assert set(calls) == set(h.LOCI)
