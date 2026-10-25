"""Dataset Hi-C sintetico con struttura nota, per test e demo senza dati reali.

Contatti = decadimento 1/s x arricchimento nei TAD x scacchiera A/B.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Union

import cooler
import numpy as np
import pandas as pd

BINSIZE = 10_000
RESOLUTIONS = [10_000, 20_000, 40_000]
CHROMSIZES = pd.Series({"chr1": 2_000_000, "chr2": 1_600_000})
TAD_BP = 200_000    # TAD da 200 kb
COMP_BP = 400_000   # blocchi A/B alternati da 400 kb, il primo è A


def _chrom_pixels(n_bins: int, offset: int) -> pd.DataFrame:
    i, j = np.triu_indices(n_bins)
    same_tad = (i * BINSIZE // TAD_BP) == (j * BINSIZE // TAD_BP)
    comp = (np.arange(n_bins) * BINSIZE // COMP_BP) % 2
    same_comp = comp[i] == comp[j]
    counts = np.round(
        2000.0 / (j - i + 1) * np.where(same_tad, 4.0, 1.0) * np.where(same_comp, 1.5, 0.6)
    ).astype(int)
    keep = counts > 0
    return pd.DataFrame({"bin1_id": i[keep] + offset, "bin2_id": j[keep] + offset,
                         "count": counts[keep]})


def make_synthetic_mcool(outdir: Union[str, Path]) -> str:
    """Scrive ``synthetic.cool`` (grezzo) e ``synthetic.mcool`` (bilanciato a 10/20/40 kb).

    Returns:
        Percorso del .mcool.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bins = cooler.binnify(CHROMSIZES, BINSIZE)

    pixels, offset = [], 0
    for chrom in CHROMSIZES.index:
        n = int((bins["chrom"] == chrom).sum())
        pixels.append(_chrom_pixels(n, offset))
        offset += n

    base = str(outdir / "synthetic.cool")
    cooler.create_cooler(base, bins, pd.concat(pixels, ignore_index=True), ordered=True)

    # scrive su un file temporaneo e rinomina solo alla fine: chi controlla
    # l'esistenza del .mcool non trova mai un file a metà (es. senza pesi)
    mcool = str(outdir / "synthetic.mcool")
    tmp = mcool + ".tmp"
    cooler.zoomify_cooler(base, tmp, RESOLUTIONS, chunksize=10_000_000)
    for res in RESOLUTIONS:
        # mad_max=0: i bin ai bordi hanno meno copertura ma non sono "cattivi".
        # Sul sintetico privo di rumore ICE non arriva a tol=1e-5: i pesi vanno bene lo stesso.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Iteration limit reached")
            cooler.balance_cooler(cooler.Cooler(f"{tmp}::/resolutions/{res}"), mad_max=0, store=True)
    Path(tmp).replace(mcool)
    return mcool


def true_compartment(starts: np.ndarray) -> np.ndarray:
    """+1 per i bin in compartimento A, -1 per B."""
    return np.where((np.asarray(starts) // COMP_BP) % 2 == 0, 1.0, -1.0)
