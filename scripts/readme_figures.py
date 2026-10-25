"""Genera le figure del README in docs/img/ da un .mcool hg38 reale.

Uso:
    uv run python scripts/readme_figures.py data/4DNFIXP4QG5B.mcool
    uv run python scripts/readme_figures.py results/synthetic/synthetic.mcool --demo
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from clew import analysis, hic_io  # noqa: E402
from clew.circos import cis_links, plot_circos, trans_links  # noqa: E402

REAL = dict(chrom="chr2", zoom="chr2:10Mb-20Mb", tad="chr2:12Mb-16Mb",
            res=dict(overview=250_000, zoom=10_000, tad=5_000, comp=100_000, ps=10_000, circos=1_000_000))
DEMO = dict(chrom="chr2", zoom="chr2:0-1000000", tad="chr2:0-800000",
            res=dict(overview=40_000, zoom=10_000, tad=10_000, comp=40_000, ps=10_000, circos=40_000))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("mcool", type=Path)
    parser.add_argument("--fasta", type=Path, default=ROOT / "data/genomes/hg38/hg38.fa")
    parser.add_argument("--out", type=Path, default=ROOT / "docs/img")
    parser.add_argument("--demo", action="store_true", help="regioni e risoluzioni del sintetico")
    args = parser.parse_args()

    cfg = DEMO if args.demo else REAL
    res, chrom, out = cfg["res"], cfg["chrom"], args.out
    out.mkdir(parents=True, exist_ok=True)
    fasta = str(args.fasta) if args.fasta.exists() and not args.demo else None
    label = args.mcool.stem

    def step(name: str) -> None:
        print(f"[{time.strftime('%H:%M:%S')}] {name}", flush=True)

    step("cromosoma intero")
    clr = hic_io.open_cooler(args.mcool, res["overview"])
    whole = (chrom, 0, int(clr.chromsizes[chrom]))
    ax = hic_io.plot_matrix(hic_io.fetch_region(args.mcool, res["overview"], whole), whole,
                            title=f"{label}: {chrom} @ {res['overview'] // 1000} kb")
    ax.figure.savefig(out / "chromosome.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    step("zoom")
    ax = hic_io.plot_matrix(hic_io.fetch_region(args.mcool, res["zoom"], cfg["zoom"]), cfg["zoom"],
                            title=f"{cfg['zoom']} @ {res['zoom'] // 1000} kb")
    ax.figure.savefig(out / "zoom.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    step("TAD e insulation")
    clr_t = hic_io.open_cooler(args.mcool, res["tad"])
    region = hic_io.cooler.util.parse_region(cfg["tad"], clr_t.chromsizes)
    ins = analysis.insulation(clr_t, regions=cfg["tad"])
    w = res["tad"] * 10
    bounds = analysis.tad_boundaries(clr_t, window_bp=w, insulation_df=ins)
    fig, (ax_m, ax_i) = plt.subplots(2, 1, figsize=(6, 7.5), sharex=True, constrained_layout=True,
                                     gridspec_kw=dict(height_ratios=[4, 1]))
    hic_io.plot_matrix(hic_io.fetch_region(args.mcool, res["tad"], region), region, ax=ax_m,
                       colorbar=False, title=f"{cfg['tad']} @ {res['tad'] // 1000} kb: TAD")
    ax_m.set_aspect("auto")
    ax_m.set_xlabel("")
    mid = (ins["start"] + ins["end"]) / 2e6
    ax_i.plot(mid, ins[f"log2_insulation_score_{w}"], lw=1, color="black")
    for b in (bounds["start"] + bounds["end"]) / 2e6:
        ax_m.axvline(b, color="tab:blue", lw=0.6, alpha=0.6)
        ax_i.axvline(b, color="tab:blue", lw=0.6, alpha=0.6)
    ax_i.set(xlabel=f"{chrom} (Mb)", ylabel=f"insulation {w // 1000} kb")
    fig.savefig(out / "tads.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    step("compartimenti")
    clr_c = hic_io.open_cooler(args.mcool, res["comp"])
    phasing = analysis.gc_phasing_track(clr_c, fasta, chrom) if fasta else None
    eigs = analysis.compartments(clr_c, phasing, chrom)
    e1 = eigs[eigs["chrom"] == chrom]
    fig, (ax_m, ax_e) = plt.subplots(2, 1, figsize=(6, 7.5), sharex=True, constrained_layout=True,
                                     gridspec_kw=dict(height_ratios=[4, 1]))
    hic_io.plot_matrix(hic_io.fetch_region(args.mcool, res["comp"], whole), whole, ax=ax_m,
                       colorbar=False, title=f"{chrom}: compartimenti A/B")
    ax_m.set_aspect("auto")
    ax_m.set_xlabel("")
    x = (e1["start"] + e1["end"]) / 2e6
    ax_e.fill_between(x, e1["E1"], where=e1["E1"] > 0, color="tab:red", step="mid")
    ax_e.fill_between(x, e1["E1"], where=e1["E1"] < 0, color="tab:blue", step="mid")
    ax_e.set(xlabel=f"{chrom} (Mb)", ylabel="E1")
    fig.savefig(out / "compartments.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    step("P(s)")
    cvd = analysis.contact_probability(hic_io.open_cooler(args.mcool, res["ps"]), chrom, nproc=4)
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    analysis.plot_contact_probability(cvd, ax=ax, label=chrom)
    fig.savefig(out / "ps.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    step("Circos")
    clr_g = hic_io.open_cooler(args.mcool, res["circos"])
    chroms = list(analysis.make_view(clr_g)["chrom"])
    phasing_g = analysis.gc_phasing_track(clr_g, fasta) if fasta else None
    links = pd.concat([trans_links(clr_g, chroms, top_n=1000),
                       cis_links(clr_g, chroms, min_distance=400_000 if args.demo else 10_000_000,
                                 top_n=300)], ignore_index=True)
    links.to_csv(ROOT / "results" / f"{label}_circos_links.tsv", sep="\t", index=False)
    plot_circos(clr_g.chromsizes[chroms], links,
                track=analysis.compartments(clr_g, phasing_g), highlight=chrom,
                out=str(out / "circos.png"), title=f"{label}: contatti trans e cis a lungo raggio")
    plt.close("all")
    step(f"fatto: {out}")


if __name__ == "__main__":
    main()
