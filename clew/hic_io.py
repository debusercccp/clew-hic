"""I/O e visualizzazione di matrici Hi-C (.cool / .mcool), sempre per regione.

Nessuna funzione carica un cromosoma intero: le matrici si leggono con
``cooler.Cooler(...).matrix().fetch(region)`` e c'è un limite esplicito
al numero di bin (``max_bins``) per evitare di saturare la RAM.

I file ``.hic`` (Juicer) non sono letti direttamente: convertili prima in
``.mcool``, per esempio con ``hictk convert in.hic out.mcool``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import cooler
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter

PathLike = Union[str, Path]
Region = Union[str, tuple[str, int, int]]

# 5000 x 5000 float64 = ~200 MB: oltre questo conviene una risoluzione più bassa.
DEFAULT_MAX_BINS = 5000

# Colormap "fall" di cooltools, copiata qui: cooltools.lib.plotting 0.7.x
# importa API rimosse in matplotlib >= 3.9 e non si può importare.
_FALL = np.array([
    (255, 255, 255), (255, 255, 204), (255, 237, 160), (254, 217, 118),
    (254, 178, 76), (253, 141, 60), (252, 78, 42), (227, 26, 28),
    (189, 0, 38), (128, 0, 38), (0, 0, 0),
]) / 255
if "fall" not in mpl.colormaps:
    mpl.colormaps.register(LinearSegmentedColormap.from_list("fall", _FALL, N=256))


def _check_not_hic(path: str) -> None:
    if path.lower().endswith(".hic"):
        raise ValueError(
            f"{path}: i file .hic non sono supportati direttamente. "
            "Convertili in .mcool, ad es. `hictk convert in.hic out.mcool`."
        )


def list_resolutions(path: PathLike, verbose: bool = True) -> list[int]:
    """Elenca le risoluzioni (bp) disponibili in un .mcool (o il binsize di un .cool).

    Args:
        path: percorso del file .mcool o .cool.
        verbose: se True stampa le risoluzioni.

    Returns:
        Lista ordinata delle risoluzioni in bp.
    """
    path = str(path)
    _check_not_hic(path)
    if cooler.fileops.is_multires_file(path):
        resolutions = sorted(
            int(group.rsplit("/", 1)[-1])
            for group in cooler.fileops.list_coolers(path)
            if group.startswith("/resolutions/")
        )
    elif cooler.fileops.is_cooler(path):
        resolutions = [cooler.Cooler(path).binsize]
    else:
        raise ValueError(f"{path} non è un file cooler valido")

    if verbose:
        print(f"{Path(path).name}: " + ", ".join(_fmt_bp(r) for r in resolutions))
    return resolutions


def open_cooler(path: PathLike, resolution: Optional[int] = None) -> cooler.Cooler:
    """Apre un Cooler alla risoluzione richiesta (senza leggere pixel).

    Args:
        path: .mcool (``resolution`` obbligatoria) o .cool.
        resolution: risoluzione in bp.

    Raises:
        ValueError: se la risoluzione non è presente nel file.
    """
    path = str(path)
    _check_not_hic(path)
    if cooler.fileops.is_multires_file(path):
        available = list_resolutions(path, verbose=False)
        if resolution not in available:
            raise ValueError(
                f"risoluzione {resolution} non disponibile in {path}; "
                f"disponibili: {available}"
            )
        return cooler.Cooler(f"{path}::/resolutions/{resolution}")

    clr = cooler.Cooler(path)
    if resolution is not None and clr.binsize != resolution:
        raise ValueError(f"{path} ha binsize {clr.binsize}, non {resolution}")
    return clr


def fetch_region(
    path: PathLike,
    resolution: int,
    region: Region,
    balance: Union[bool, str] = True,
    max_bins: int = DEFAULT_MAX_BINS,
) -> np.ndarray:
    """Ritorna la matrice di contatto (densa) di una regione cis.

    Args:
        path: file .mcool o .cool.
        resolution: risoluzione in bp.
        region: es. ``"chr2:10,000,000-20,000,000"``, ``"chr2"`` o
            ``("chr2", 10_000_000, 20_000_000)``.
        balance: True usa la colonna ``weight``; una stringa sceglie
            un'altra colonna di pesi; False ritorna i conteggi grezzi.
        max_bins: numero massimo di bin per lato prima di rifiutare la lettura.

    Returns:
        Matrice quadrata ``float`` (NaN sui bin filtrati dal bilanciamento).
    """
    clr = open_cooler(path, resolution)
    chrom, start, end = cooler.util.parse_region(region, clr.chromsizes)

    lo, hi = clr.extent((chrom, start, end))
    if hi - lo > max_bins:
        raise ValueError(
            f"{chrom}:{start}-{end} a {_fmt_bp(resolution)} sono {hi - lo} bin "
            f"(limite {max_bins}): restringi la regione o usa una risoluzione più bassa"
        )

    weight_name = "weight" if balance is True else balance
    if weight_name and weight_name not in clr.bins().columns:
        raise ValueError(
            f"la colonna di pesi '{weight_name}' manca: il file non è bilanciato "
            "(usa `cooler balance` oppure balance=False)"
        )
    return clr.matrix(balance=balance).fetch((chrom, start, end))


def plot_matrix(
    mat: np.ndarray,
    region: Region,
    log: bool = True,
    cmap: str = "fall",
    out: Optional[PathLike] = None,
    ax: Optional[Axes] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    colorbar: bool = True,
    title: Optional[str] = None,
) -> Axes:
    """Disegna una heatmap di contatti con assi in Mb.

    Args:
        mat: matrice quadrata, tipicamente da :func:`fetch_region`.
        region: la regione usata per ``mat`` (servono coordinate esplicite).
        log: scala colori logaritmica (zeri e NaN vengono mascherati).
        cmap: nome della colormap; "fall" è quella di cooltools.
        out: se dato, salva la figura (png, pdf, svg...).
        ax: Axes su cui disegnare; se None ne crea uno nuovo.
        vmin, vmax: limiti della scala colori (default: min/max dei dati validi).
        colorbar: aggiunge la colorbar.
        title: titolo; default la regione.

    Returns:
        L'Axes con la heatmap.
    """
    chrom, start, end = cooler.util.parse_region(region)
    data = np.ma.masked_invalid(np.asarray(mat, dtype=float))
    if log:
        data = np.ma.masked_less_equal(data, 0)

    valid = data.compressed()
    if valid.size:
        vmin = valid.min() if vmin is None else vmin
        vmax = valid.max() if vmax is None else vmax
    norm = LogNorm(vmin, vmax) if log and valid.size else Normalize(vmin, vmax)

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5.5), constrained_layout=True)

    # L'estensione usa le coordinate della regione: l'errore è al più un bin.
    extent = (start / 1e6, end / 1e6, end / 1e6, start / 1e6)
    im = ax.imshow(data, cmap=cmap, norm=norm, extent=extent, interpolation="none")

    mb = FuncFormatter(lambda x, _: f"{x:g}")
    ax.xaxis.set_major_formatter(mb)
    ax.yaxis.set_major_formatter(mb)
    ax.set_xlabel(f"{chrom} (Mb)")
    ax.set_ylabel(f"{chrom} (Mb)")
    ax.set_title(title or f"{chrom}:{start:,}-{end:,}")

    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                           label="contatti" + (" (log)" if log else ""))
    if out is not None:
        ax.figure.savefig(out, dpi=150, bbox_inches="tight")
    return ax


def plot_pair(
    path: PathLike,
    resolution: int,
    region_a: Region,
    region_b: Region,
    balance: Union[bool, str] = True,
    log: bool = True,
    cmap: str = "fall",
    out: Optional[PathLike] = None,
) -> Figure:
    """Confronta due regioni dello stesso file affiancate, con scala colori comune.

    Args:
        path: file .mcool o .cool.
        resolution: risoluzione in bp.
        region_a, region_b: regioni cis da confrontare.
        balance, log, cmap: come in :func:`fetch_region` / :func:`plot_matrix`.
        out: se dato, salva la figura.

    Returns:
        La Figure con i due pannelli e una colorbar condivisa.
    """
    mats = [fetch_region(path, resolution, r, balance=balance) for r in (region_a, region_b)]

    valid = np.concatenate([m[np.isfinite(m) & ((m > 0) if log else True)] for m in mats])
    vmin, vmax = (valid.min(), valid.max()) if valid.size else (None, None)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True)
    for ax, mat, region in zip(axes, mats, (region_a, region_b)):
        plot_matrix(mat, region, log=log, cmap=cmap, ax=ax,
                    vmin=vmin, vmax=vmax, colorbar=False)
    fig.colorbar(axes[1].images[0], ax=axes, shrink=0.8,
                 label="contatti" + (" (log)" if log else ""))
    fig.suptitle(f"{Path(str(path)).name} @ {_fmt_bp(resolution)}")

    if out is not None:
        fig.savefig(out, dpi=150, bbox_inches="tight")
    return fig


def _fmt_bp(bp: int) -> str:
    if bp >= 1_000_000 and bp % 100_000 == 0:
        return f"{bp / 1_000_000:g} Mb"
    if bp >= 1_000 and bp % 1_000 == 0:
        return f"{bp // 1_000} kb"
    return f"{bp} bp"
