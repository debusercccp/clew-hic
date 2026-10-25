import numpy as np
import pandas as pd
import pytest

from clew import analysis, hic_io
from tests.conftest import CHROMSIZES, COMP_BP, TAD_BP, true_compartment


@pytest.fixture(scope="module")
def clr10(synthetic_mcool):
    return hic_io.open_cooler(synthetic_mcool, 10_000)


def test_make_view_regions(clr10):
    view = analysis.make_view(clr10, ["chr1", "chr2:100,000-500,000"])
    assert view[["chrom", "start", "end"]].values.tolist() == [
        ["chr1", 0, 2_000_000], ["chr2", 100_000, 500_000]]


def test_contact_probability_decays(clr10):
    cvd = analysis.contact_probability(clr10, ["chr1"])
    assert isinstance(cvd, pd.DataFrame)
    assert cvd["s_bp"].min() == 20_000
    p = cvd.dropna(subset=["P"])
    assert np.corrcoef(np.log(p["s_bp"]), np.log(p["P"]))[0, 1] < -0.9
    # oltre ~1/4 del cromosoma sintetico restano pochi pixel e domina la scacchiera A/B
    short = p[p["s_bp"] <= 500_000]
    assert (short["slope"] < 0).all()


def test_compartments_recover_ab(synthetic_mcool):
    clr = hic_io.open_cooler(synthetic_mcool, 40_000)
    bins = clr.bins()[:][["chrom", "start", "end"]]
    phasing = bins.assign(value=true_compartment(bins["start"]))

    eigs = analysis.compartments(clr, phasing_track=phasing)
    ok = eigs["E1"].notna()
    assert ok.mean() > 0.9
    truth = np.where(true_compartment(eigs["start"]) > 0, "A", "B")
    assert (eigs.loc[ok, "compartment"] == truth[ok]).mean() > 0.95
    assert len(eigs.attrs["eigvals"]) == len(CHROMSIZES)


def test_compartments_warns_without_phasing(synthetic_mcool):
    clr = hic_io.open_cooler(synthetic_mcool, 40_000)
    with pytest.warns(UserWarning, match="arbitrario"):
        analysis.compartments(clr, regions="chr1")


def test_insulation_rejects_bad_window(clr10):
    with pytest.raises(ValueError, match="multipli"):
        analysis.insulation(clr10, [15_000])


def test_tad_boundaries_match_synthetic_tads(clr10):
    # threshold=0: tutti i minimi locali. Con "Li" restano solo i confini più
    # forti (quelli che sono anche cambi di compartimento), ed è corretto così.
    bounds = analysis.tad_boundaries(clr10, window_bp=100_000, threshold=0)
    assert {"chrom", "start", "insulation", "strength"} <= set(bounds.columns)

    # confini veri: multipli di 200 kb interni a ogni cromosoma
    truth = {(c, s) for c, size in CHROMSIZES.items() for s in range(TAD_BP, size, TAD_BP)}
    near_truth = [
        any(c == chrom and abs(s - start) <= 10_000 for c, s in truth)
        for chrom, start in zip(bounds["chrom"], bounds["start"])
    ]
    assert all(near_truth), bounds[~np.array(near_truth)]
    assert len(bounds) >= 0.8 * len(truth)


def test_tad_boundaries_li_keeps_compartment_switches(clr10):
    bounds = analysis.tad_boundaries(clr10, window_bp=100_000)
    switches = {(c, s) for c, size in CHROMSIZES.items() for s in range(COMP_BP, size, COMP_BP)}
    assert len(bounds) > 0
    for chrom, start in zip(bounds["chrom"], bounds["start"]):
        assert any(c == chrom and abs(s - start) <= 10_000 for c, s in switches)


def test_cli_writes_tsv(synthetic_mcool, tmp_path):
    out = tmp_path / "bounds.tsv"
    analysis._main(["boundaries", synthetic_mcool, "--res", "10000",
                    "--region", "chr1", "--window", "100000", "--out", str(out)])
    df = pd.read_csv(out, sep="\t")
    assert len(df) > 0 and set(df["chrom"]) == {"chr1"}


def test_gc_phasing_track_reads_fasta(synthetic_mcool, tmp_path):
    # genoma finto: blocchi A (GC) e B (AT) alternati come i compartimenti sintetici
    fasta = tmp_path / "toy.fa"
    with open(fasta, "w") as fh:
        for chrom, size in CHROMSIZES.items():
            seq = "".join(("GC" if (i // COMP_BP) % 2 == 0 else "AT") * (COMP_BP // 2)
                          for i in range(0, size, COMP_BP))[:size]
            fh.write(f">{chrom}\n" + "\n".join(seq[k:k + 80] for k in range(0, size, 80)) + "\n")

    clr = hic_io.open_cooler(synthetic_mcool, 40_000)
    gc = analysis.gc_phasing_track(clr, str(fasta), "chr2")
    assert set(gc["chrom"]) == {"chr2"}
    assert np.allclose(gc["GC"], (true_compartment(gc["start"]) > 0).astype(float))
