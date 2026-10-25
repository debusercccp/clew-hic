"""Vista circolare (stile Circos) dei contatti tra cromosomi.

I cromosomi sono archi su un cerchio; le corde collegano coppie di bin con
contatti più alti dell'atteso: trans (cromosomi diversi, magenta) e, se
richieste, cis a lungo raggio (stesso cromosoma, arancione). Un anello
opzionale mostra una traccia per bin (es. E1 dei compartimenti: A rosso, B blu).

Serve il genoma intero, quindi solo a bassa risoluzione (>= 500 kb):
a 1 Mb la matrice hg38 ha ~3 100 bin (~76 MB).
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Optional, Union

import cooler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MplPath

from .analysis import make_view
from .hic_io import DEFAULT_MAX_BINS

GAP_DEG = 1.5  # spazio tra un cromosoma e l'altro


def _usable_bins(
    clr: cooler.Cooler,
    chroms: Optional[Sequence[str]],
    edge_bins: int,
    min_valid_frac: float,
    max_bins: int,
) -> dict[str, np.ndarray]:
    """Maschera per cromosoma dei bin affidabili per le corde.

    Le corde più "arricchite" cadono quasi sempre su regioni ripetute dove le
    letture mappano male: bordi dei centromeri e dei buchi dell'assembly,
    estremità dei cromosomi, cromosomi quasi vuoti (chrY in una linea
    femminile). Si escludono i cromosomi con meno di ``min_valid_frac`` bin
    validi e i bin entro ``edge_bins`` da un bin filtrato o da un'estremità.
    """
    chroms = list(chroms) if chroms is not None else list(make_view(clr)["chrom"])
    n = sum(hi - lo for lo, hi in (clr.extent(c) for c in chroms))
    if n > max_bins:
        raise ValueError(f"{n} bin (limite {max_bins}): usa una risoluzione più bassa (>= 500 kb)")

    usable = {}
    for chrom in chroms:
        valid = clr.bins().fetch(chrom)["weight"].notna().to_numpy()
        if valid.mean() < min_valid_frac:
            continue
        bad = np.pad(~valid, edge_bins, constant_values=True)  # le estremità contano come buchi
        near_bad = np.convolve(bad, np.ones(2 * edge_bins + 1), mode="valid") > 0
        usable[chrom] = valid & ~near_bad
    return usable


def trans_links(
    clr: cooler.Cooler,
    chroms: Optional[Sequence[str]] = None,
    top_n: int = 200,
    edge_bins: int = 2,
    min_valid_frac: float = 0.5,
    max_bins: int = DEFAULT_MAX_BINS,
) -> pd.DataFrame:
    """Coppie di bin trans più arricchite rispetto alla media della coppia di cromosomi.

    Args:
        clr: Cooler bilanciato a bassa risoluzione (>= 500 kb).
        chroms: cromosomi da includere (default: tutti tranne mitocondrio e contig).
        top_n: quante coppie tenere, dalla più arricchita.
        edge_bins: esclude i bin a questa distanza da centromeri, buchi ed estremità.
        min_valid_frac: esclude i cromosomi con meno di questa frazione di bin validi.
        max_bins: limite sui bin totali della matrice letta.

    Returns:
        chrom1, start1, end1, chrom2, start2, end2, balanced, oe (observed/expected).
    """
    usable = _usable_bins(clr, chroms, edge_bins, min_valid_frac, max_bins)
    chroms = list(usable)
    extents = {c: clr.extent(c) for c in chroms}

    bins = clr.bins()[:]
    rows = []
    for a, ca in enumerate(chroms):
        for cb in chroms[a + 1:]:
            block = clr.matrix(balance=True).fetch(ca, cb)
            block[~usable[ca], :] = np.nan
            block[:, ~usable[cb]] = np.nan
            expected = np.nanmean(block)
            if not np.isfinite(expected) or expected <= 0:
                continue
            i, j = np.nonzero(np.isfinite(block) & (block > 0))
            rows.append(pd.DataFrame({
                "bin1": i + extents[ca][0], "bin2": j + extents[cb][0],
                "balanced": block[i, j], "oe": block[i, j] / expected,
            }))

    return _top_links(rows, bins, top_n)


def cis_links(
    clr: cooler.Cooler,
    chroms: Optional[Sequence[str]] = None,
    min_distance: int = 10_000_000,
    top_n: int = 200,
    edge_bins: int = 2,
    min_valid_frac: float = 0.5,
    max_bins: int = DEFAULT_MAX_BINS,
) -> pd.DataFrame:
    """Coppie di bin cis a lungo raggio più arricchite rispetto all'atteso per la loro distanza.

    L'atteso di ogni cella è la media della sua diagonale (stessa distanza
    genomica, stesso cromosoma), quindi il decadimento con la distanza non conta.

    Args:
        clr: Cooler bilanciato a bassa risoluzione (>= 500 kb).
        chroms: cromosomi da includere (default: tutti tranne mitocondrio e contig).
        min_distance: distanza minima in bp tra i due bin.
        top_n: quante coppie tenere, dalla più arricchita.
        edge_bins, min_valid_frac: filtri sulle regioni inaffidabili, come in :func:`trans_links`.
        max_bins: limite sui bin totali.

    Returns:
        Stesse colonne di :func:`trans_links`, con chrom1 == chrom2.
    """
    usable = _usable_bins(clr, chroms, edge_bins, min_valid_frac, max_bins)

    min_diag = max(1, -(-min_distance // clr.binsize))
    bins = clr.bins()[:]
    rows = []
    for chrom, ok_bins in usable.items():
        mat = clr.matrix(balance=True).fetch(chrom)
        mat[~ok_bins, :] = np.nan
        mat[:, ~ok_bins] = np.nan
        lo = clr.extent(chrom)[0]
        for d in range(min_diag, mat.shape[0]):
            diag = np.diagonal(mat, d)
            ok = np.isfinite(diag) & (diag > 0)
            expected = diag[ok].mean() if ok.any() else np.nan
            if not np.isfinite(expected) or expected <= 0:
                continue
            i = np.flatnonzero(ok)
            rows.append(pd.DataFrame({
                "bin1": i + lo, "bin2": i + d + lo,
                "balanced": diag[i], "oe": diag[i] / expected,
            }))
    return _top_links(rows, bins, top_n)


_LINK_COLUMNS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "balanced", "oe"]


def _top_links(rows: list[pd.DataFrame], bins: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Tiene le ``top_n`` coppie per oe e aggiunge le coordinate dei bin."""
    if not rows:
        return pd.DataFrame(columns=_LINK_COLUMNS)
    links = pd.concat(rows, ignore_index=True).nlargest(top_n, "oe")
    b1 = bins.loc[links["bin1"], ["chrom", "start", "end"]].to_numpy()
    b2 = bins.loc[links["bin2"], ["chrom", "start", "end"]].to_numpy()
    out = pd.DataFrame(np.hstack([b1, b2]), columns=_LINK_COLUMNS[:6])
    out["chrom1"] = out["chrom1"].astype(str)
    out["chrom2"] = out["chrom2"].astype(str)
    out[["start1", "end1", "start2", "end2"]] = out[["start1", "end1", "start2", "end2"]].astype(int)
    out["balanced"] = links["balanced"].to_numpy()
    out["oe"] = links["oe"].to_numpy()
    return out.reset_index(drop=True)


class _Layout:
    """Converte (cromosoma, posizione) in angolo, con i cromosomi in senso orario da ore 12."""

    def __init__(self, chromsizes: pd.Series):
        self.sizes = chromsizes
        span = 360.0 - GAP_DEG * len(chromsizes)
        self.deg_per_bp = span / chromsizes.sum()
        starts = np.concatenate([[0], np.cumsum(chromsizes.to_numpy() * self.deg_per_bp + GAP_DEG)[:-1]])
        self.start = dict(zip(chromsizes.index, starts))

    def angle(self, chrom: str, pos: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        # 90° = ore 12, poi in senso orario
        return 90.0 - (self.start[chrom] + np.asarray(pos) * self.deg_per_bp)


def _xy(deg: float, r: float) -> tuple[float, float]:
    t = np.deg2rad(deg)
    return r * np.cos(t), r * np.sin(t)


def plot_circos(
    chromsizes: pd.Series,
    links: Optional[pd.DataFrame] = None,
    track: Optional[pd.DataFrame] = None,
    track_col: str = "E1",
    highlight: Optional[str] = None,
    link_color: str = "#b5179e",
    cis_color: str = "#f28e2b",
    highlight_color: str = "#2ca02c",
    ax: Optional[Axes] = None,
    out: Optional[str] = None,
    title: Optional[str] = None,
) -> Axes:
    """Disegna il grafico circolare.

    Args:
        chromsizes: dimensioni dei cromosomi da mostrare, nell'ordine voluto.
        links: output di :func:`trans_links` e/o :func:`cis_links` (anche concatenati);
            spessore e opacità crescono con ``oe``, scalato separatamente per cis e trans.
        track: DataFrame per bin (chrom, start, end, ``track_col``), es. da
            ``analysis.compartments``; valori > 0 in rosso, < 0 in blu.
        track_col: colonna di ``track`` da disegnare.
        highlight: regione (es. ``"chr2"`` o ``"chr2:10Mb-20Mb"``): le corde che
            la toccano sono disegnate sopra le altre, in ``highlight_color``.
        link_color: colore delle corde trans (stile Circos classico: un solo colore).
        cis_color: colore delle corde cis a lungo raggio.
        ax: Axes su cui disegnare; se None ne crea uno.
        out: se dato, salva la figura.
        title: titolo della figura.

    Returns:
        L'Axes con il grafico.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    lay = _Layout(chromsizes)

    r_ideo, w_ideo = 1.0, 0.06            # anello dei cromosomi
    r_track, w_track = 0.92, 0.12         # anello della traccia (barre attorno a 0.86)
    r_link = 0.79 if track is not None else 0.93

    # anello dei cromosomi: segmenti alternati bianco/nero con bordo, come in Circos
    for k, (chrom, size) in enumerate(chromsizes.items()):
        a0, a1 = lay.angle(chrom, size), lay.angle(chrom, 0)
        ax.add_patch(Wedge((0, 0), r_ideo, a0, a1, width=w_ideo,
                           fc="black" if k % 2 == 0 else "white", ec="black", lw=0.8))
        mid = (a0 + a1) / 2
        x, y = _xy(mid, r_ideo + 0.06)
        ax.text(x, y, chrom.removeprefix("chr"), ha="center", va="center",
                fontsize=8, rotation=mid - 90 if 0 <= (mid % 360) <= 180 else mid + 90)

    if track is not None:
        values = track[track["chrom"].isin(chromsizes.index)].dropna(subset=[track_col])
        vmax = np.nanmax(np.abs(values[track_col])) or 1.0
        base = r_track - w_track / 2
        for chrom, start, end, v in values[["chrom", "start", "end", track_col]].itertuples(index=False):
            h = (w_track / 2) * v / vmax
            ax.add_patch(Wedge((0, 0), base + max(h, 0), lay.angle(chrom, end), lay.angle(chrom, start),
                               width=abs(h), color="tab:red" if v > 0 else "tab:blue", lw=0))

    if links is not None and len(links):
        links = links[links["chrom1"].isin(chromsizes.index) & links["chrom2"].isin(chromsizes.index)]
        hit = np.zeros(len(links), dtype=bool)
        if highlight:
            hc, hs, he = cooler.util.parse_region(highlight, chromsizes)
            for c, s0, s1 in (("chrom1", "start1", "end1"), ("chrom2", "start2", "end2")):
                hit |= ((links[c] == hc) & (links[s1] > hs) & (links[s0] < he)).to_numpy()
        cis = (links["chrom1"] == links["chrom2"]).to_numpy()
        # oe di cis e trans hanno scale diverse: si normalizzano separatamente
        strength = np.zeros(len(links))
        for mask in (cis, ~cis):
            if mask.any():
                oe = links["oe"].to_numpy(dtype=float)[mask]
                strength[mask] = (oe - oe.min()) / (np.ptp(oe) or 1.0)
        # prima le corde normali, poi quelle evidenziate sopra
        for highlighted in (False, True):
            sel = hit == highlighted
            for row, s, is_cis in zip(links[sel].itertuples(index=False), strength[sel], cis[sel]):
                a1 = lay.angle(row.chrom1, (row.start1 + row.end1) / 2)
                a2 = lay.angle(row.chrom2, (row.start2 + row.end2) / 2)
                # trans: curva che passa vicino al centro; cis: "petalo" verso l'interno,
                # più profondo quanto più i due punti sono distanti
                depth = max(0.35, 0.7 - abs(a1 - a2) / 90)
                ctrl = (0.0, 0.0) if not is_cis else _xy((a1 + a2) / 2, r_link * depth)
                path = MplPath([_xy(a1, r_link), ctrl, _xy(a2, r_link)],
                               [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
                ax.add_patch(PathPatch(
                    path, fc="none",
                    ec=highlight_color if highlighted else (cis_color if is_cis else link_color),
                    lw=1.1 if highlighted else (0.8 + 0.8 * s if is_cis else 0.5 + 0.8 * s),
                    alpha=0.95 if highlighted else (0.5 + 0.4 * s if is_cis else 0.2 + 0.4 * s),
                    zorder=3.5 if highlighted else (2.5 if is_cis else 2)))

    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title)
    if out is not None:
        ax.figure.savefig(out, dpi=150, bbox_inches="tight")
    return ax
