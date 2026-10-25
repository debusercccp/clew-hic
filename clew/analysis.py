"""Analisi Hi-C con cooltools: P(s), compartimenti A/B, insulation e confini di TAD.

Ogni funzione accetta un ``cooler.Cooler`` (vedi ``hic_io.open_cooler``) e
ritorna un ``pandas.DataFrame``. Il parametro ``regions`` limita l'analisi a
cromosomi o sotto-regioni (es. ``["chr2"]`` o ``["chr2:10Mb-60Mb"]``).

Uso da riga di comando (comodo nel servizio ``worker`` di docker compose)::

    python -m clew.analysis insulation data/x.mcool --res 10000 --region chr2 \\
        --out results/x_insulation.tsv
"""
from __future__ import annotations

import argparse
import warnings
from collections.abc import Sequence
from typing import Optional, Union

import bioframe
import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

Regions = Optional[Union[str, Sequence[str]]]

_SKIP_CHROMS = {"chrM", "M", "chrMT", "MT"}

# Pixel letti per blocco da ogni processo. I default di cooltools (10-20 M)
# costano ~4,6 GB per processo su un file da 4 miliardi di contatti; con 2 M
# si scende a ~1,1 GB (misurato su GM12878 4DNFIXP4QG5B, chr2 @ 10 kb) per
# un tempo appena maggiore (+11%).
CHUNKSIZE = 2_000_000


def make_view(clr: cooler.Cooler, regions: Regions = None) -> pd.DataFrame:
    """Costruisce la view (chrom, start, end, name) richiesta da cooltools.

    Args:
        clr: Cooler di riferimento (per le dimensioni dei cromosomi).
        regions: None = tutti i cromosomi tranne mitocondrio e contig
            non posizionati; altrimenti nomi di cromosoma o regioni UCSC.
    """
    if regions is None:
        regions = [c for c in clr.chromnames if c not in _SKIP_CHROMS and "_" not in c]
    elif isinstance(regions, str):
        regions = [regions]

    rows = []
    for region in regions:
        chrom, start, end = cooler.util.parse_region(region, clr.chromsizes)
        rows.append((chrom, start, end, f"{chrom}:{start}-{end}"))
    return bioframe.make_viewframe(pd.DataFrame(rows, columns=["chrom", "start", "end", "name"]))


def contact_probability(
    clr: cooler.Cooler,
    regions: Regions = None,
    ignore_diags: int = 2,
    nproc: int = 1,
    chunksize: int = CHUNKSIZE,
) -> pd.DataFrame:
    """Curva P(s): probabilità di contatto media in funzione della distanza.

    La RAM di picco è circa ``nproc`` x 1,1 GB con il ``chunksize`` di default.

    Returns:
        Output di ``cooltools.expected_cis`` più le colonne ``s_bp``
        (distanza in bp), ``P`` (curva lisciata e aggregata sulle regioni)
        e ``slope`` (d log P / d log s). Le diagonali ignorate sono escluse.
    """
    cvd = cooltools.expected_cis(
        clr,
        view_df=make_view(clr, regions),
        ignore_diags=ignore_diags,
        smooth=True,
        aggregate_smoothed=True,
        nproc=nproc,
        chunksize=chunksize,
    )
    cvd = cvd[cvd["dist"] >= ignore_diags].copy()
    cvd["s_bp"] = cvd["dist"] * clr.binsize
    cvd["P"] = cvd["balanced.avg.smoothed.agg"]

    cvd["slope"] = np.nan
    for _, group in cvd.groupby("region1", observed=True):
        ok = group.index[group["P"] > 0]
        if len(ok) > 1:
            cvd.loc[ok, "slope"] = np.gradient(
                np.log(cvd.loc[ok, "P"]), np.log(cvd.loc[ok, "s_bp"])
            )
    return cvd.reset_index(drop=True)


def plot_contact_probability(cvd: pd.DataFrame, ax: Optional[Axes] = None, label: Optional[str] = None) -> Axes:
    """Disegna P(s) in log-log da :func:`contact_probability` (una curva per regione)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    for region, group in cvd.groupby("region1"):
        ax.loglog(group["s_bp"], group["P"], label=label or region)
    ax.set_xlabel("distanza genomica s (bp)")
    ax.set_ylabel("P(s)")
    ax.legend(fontsize="small")
    return ax


def gc_phasing_track(clr: cooler.Cooler, fasta: str, regions: Regions = None) -> pd.DataFrame:
    """Contenuto GC per bin, da usare come ``phasing_track`` in :func:`compartments`.

    Args:
        clr: Cooler alla risoluzione dei compartimenti.
        fasta: genoma di riferimento (es. hg38.fa) dello stesso assembly del file.
        regions: limita il calcolo ai cromosomi di queste regioni (più veloce).
    """
    # bioframe >= 0.8 usa pysam di default; pyfaidx è già tra le dipendenze
    genome = bioframe.load_fasta(fasta, engine="pyfaidx")
    chroms = set(make_view(clr, regions)["chrom"]) & set(genome.keys())
    bins = clr.bins()[:][["chrom", "start", "end"]]
    bins = bins[bins["chrom"].isin(chroms)]
    return bioframe.frac_gc(bins, genome)


def compartments(
    clr: cooler.Cooler,
    phasing_track: Optional[pd.DataFrame] = None,
    regions: Regions = None,
    n_eigs: int = 3,
) -> pd.DataFrame:
    """Compartimenti A/B dal primo autovettore della matrice O/E cis.

    Nota: ``eigs_cis`` legge per intero la matrice di ogni regione della view,
    quindi va usata a bassa risoluzione (100 kb-250 kb), oppure su sotto-regioni.

    Args:
        clr: Cooler bilanciato.
        phasing_track: DataFrame (chrom, start, end, valore), es. da
            :func:`gc_phasing_track`; orienta E1 in modo che E1 > 0 = A
            (ricco in GC). Senza, il segno di E1 è arbitrario.
        regions: cromosomi/regioni da analizzare.
        n_eigs: numero di autovettori da calcolare.

    Returns:
        Bin con colonne E1..En e ``compartment`` ("A", "B" o None se NaN).
        Gli autovalori sono in ``df.attrs["eigvals"]``.
    """
    if clr.binsize < 50_000:
        warnings.warn(
            f"binsize {clr.binsize}: i compartimenti si calcolano di solito a 100-250 kb",
            stacklevel=2,
        )
    if phasing_track is None:
        warnings.warn("nessun phasing_track: il segno di E1 (A/B) è arbitrario", stacklevel=2)

    eigvals, eigvecs = cooltools.eigs_cis(
        clr, phasing_track=phasing_track, view_df=make_view(clr, regions), n_eigs=n_eigs
    )
    eigvecs["compartment"] = np.where(
        eigvecs["E1"] > 0, "A", np.where(eigvecs["E1"] < 0, "B", None)
    )
    eigvecs.attrs["eigvals"] = eigvals
    return eigvecs


def insulation(
    clr: cooler.Cooler,
    window_bp: Optional[Sequence[int]] = None,
    regions: Regions = None,
    threshold: Union[str, float] = "Li",
    nproc: int = 1,
    chunksize: int = CHUNKSIZE,
) -> pd.DataFrame:
    """Insulation score a diamante (cooltools), calcolato a blocchi.

    Args:
        clr: Cooler bilanciato (consigliati 5-25 kb).
        window_bp: finestre in bp, multipli del binsize (default 5, 10, 25 bin).
        regions: cromosomi/regioni da analizzare.
        threshold: soglia di forza dei confini ("Li", "Otsu" o un numero).

    Returns:
        Per ogni bin e finestra ``w``: ``log2_insulation_score_w``,
        ``boundary_strength_w``, ``is_boundary_w``.
    """
    binsize = clr.binsize
    if window_bp is None:
        window_bp = [binsize * k for k in (5, 10, 25)]
    window_bp = [int(w) for w in window_bp]
    bad = [w for w in window_bp if w % binsize]
    if bad:
        raise ValueError(f"le finestre {bad} non sono multipli del binsize {binsize}")

    return cooltools.insulation(
        clr, window_bp, view_df=make_view(clr, regions), threshold=threshold, nproc=nproc,
        chunksize=chunksize,
    )


def tad_boundaries(
    clr: cooler.Cooler,
    window_bp: Optional[int] = None,
    regions: Regions = None,
    min_strength: Optional[float] = None,
    threshold: Union[str, float] = "Li",
    insulation_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Confini di TAD: minimi locali dell'insulation sopra la soglia di forza.

    Args:
        clr: Cooler bilanciato.
        window_bp: finestra in bp (default 10 bin).
        regions: cromosomi/regioni (ignorato se si passa ``insulation_df``).
        min_strength: filtro aggiuntivo su ``boundary_strength``.
        threshold: soglia automatica ("Li", "Otsu") o numerica per
            ``is_boundary``; 0 tiene tutti i minimi locali.
        insulation_df: output di :func:`insulation` già calcolato, per non rifarlo.

    Returns:
        Un confine per riga: chrom, start, end, region, insulation, strength.
    """
    window_bp = int(window_bp or clr.binsize * 10)
    if insulation_df is None:
        insulation_df = insulation(clr, [window_bp], regions=regions, threshold=threshold)

    flag = f"is_boundary_{window_bp}"
    if flag not in insulation_df:
        raise ValueError(f"insulation_df non contiene la finestra {window_bp}")

    out = insulation_df.loc[
        insulation_df[flag],
        ["chrom", "start", "end", "region",
         f"log2_insulation_score_{window_bp}", f"boundary_strength_{window_bp}"],
    ].rename(columns={
        f"log2_insulation_score_{window_bp}": "insulation",
        f"boundary_strength_{window_bp}": "strength",
    })
    if min_strength is not None:
        out = out[out["strength"] >= min_strength]
    return out.reset_index(drop=True)


def _main(argv: Optional[Sequence[str]] = None) -> None:
    from .hic_io import open_cooler

    parser = argparse.ArgumentParser(prog="python -m clew.analysis", description=__doc__.split("\n")[0])
    parser.add_argument("analysis", choices=["ps", "compartments", "insulation", "boundaries"])
    parser.add_argument("path", help="file .mcool o .cool")
    parser.add_argument("--res", type=int, required=True, help="risoluzione in bp")
    parser.add_argument("--region", action="append", help="cromosoma o regione (ripetibile)")
    parser.add_argument("--window", type=int, action="append", help="finestra insulation in bp")
    parser.add_argument("--fasta", help="genoma per orientare i compartimenti con il GC")
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--out", help="file TSV di output (default: stampa le prime righe)")
    args = parser.parse_args(argv)

    clr = open_cooler(args.path, args.res)
    if args.analysis == "ps":
        df = contact_probability(clr, args.region, nproc=args.nproc)
    elif args.analysis == "compartments":
        phasing = gc_phasing_track(clr, args.fasta, args.region) if args.fasta else None
        df = compartments(clr, phasing, args.region)
    elif args.analysis == "insulation":
        df = insulation(clr, args.window, args.region, nproc=args.nproc)
    else:
        df = tad_boundaries(clr, args.window[0] if args.window else None, args.region)

    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"{len(df)} righe -> {args.out}")
    else:
        print(df.head(20).to_string())


if __name__ == "__main__":
    _main()
