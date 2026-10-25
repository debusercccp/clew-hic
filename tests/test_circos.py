import matplotlib
import matplotlib.colors
import matplotlib.patches

matplotlib.use("Agg")

import cooler
import numpy as np
import pandas as pd
import pytest

from clew.circos import _Layout, cis_links, plot_circos, trans_links

BINSIZE = 1_000_000
CHROMSIZES = pd.Series({"chr1": 20_000_000, "chr2": 15_000_000, "chr3": 10_000_000})
HOTSPOT = (("chr1", 5), ("chr3", 7))  # coppia di bin trans con contatti arricchiti
CIS_HOTSPOT = ("chr1", 2, 15)         # bin 2 e 15 di chr1: contatto cis a 13 Mb arricchito


@pytest.fixture(scope="module")
def trans_cooler(tmp_path_factory) -> cooler.Cooler:
    """Cooler a 1 Mb: cis con decadimento, trans uniforme, un hotspot trans noto."""
    bins = cooler.binnify(CHROMSIZES, BINSIZE)
    n = len(bins)
    chrom = bins["chrom"].to_numpy()
    offset = {c: int(np.flatnonzero(chrom == c)[0]) for c in CHROMSIZES.index}
    i, j = np.triu_indices(n)
    counts = np.where(chrom[i] == chrom[j], 1000 // (j - i + 1), 20)
    a = offset[HOTSPOT[0][0]] + HOTSPOT[0][1]
    b = offset[HOTSPOT[1][0]] + HOTSPOT[1][1]
    counts[(i == a) & (j == b)] = 400
    c, x, y = CIS_HOTSPOT
    counts[(i == offset[c] + x) & (j == offset[c] + y)] = 300
    pixels = pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": counts})

    path = str(tmp_path_factory.mktemp("circos") / "trans.cool")
    cooler.create_cooler(path, bins, pixels, ordered=True)
    clr = cooler.Cooler(path)
    cooler.balance_cooler(clr, mad_max=0, store=True)
    return clr


def test_trans_links_finds_hotspot(trans_cooler):
    links = trans_links(trans_cooler, top_n=5)
    assert list(links.columns) == ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "balanced", "oe"]
    top = links.iloc[0]
    assert (top.chrom1, top.start1) == (HOTSPOT[0][0], HOTSPOT[0][1] * BINSIZE)
    assert (top.chrom2, top.start2) == (HOTSPOT[1][0], HOTSPOT[1][1] * BINSIZE)
    assert top.oe > 5 and links["oe"].is_monotonic_decreasing
    assert (links["chrom1"] != links["chrom2"]).all()  # solo trans


def test_trans_links_guard(trans_cooler):
    with pytest.raises(ValueError, match="risoluzione più bassa"):
        trans_links(trans_cooler, max_bins=10)


def test_trans_links_no_trans_contacts(synthetic_mcool):
    clr = cooler.Cooler(f"{synthetic_mcool}::/resolutions/40000")  # il sintetico è solo cis
    assert trans_links(clr).empty


def test_layout_is_clockwise_from_noon():
    lay = _Layout(CHROMSIZES)
    assert lay.angle("chr1", 0) == pytest.approx(90.0)
    assert lay.angle("chr1", 10_000_000) < 90.0               # orario: l'angolo scende
    assert lay.angle("chr2", 0) < lay.angle("chr1", 20_000_000)  # c'è uno spazio tra cromosomi
    end = lay.angle("chr3", 10_000_000)
    assert end == pytest.approx(90.0 - 360.0 + 1.5)           # torna quasi a ore 12


def test_plot_circos_draws_ring_track_and_links(trans_cooler, tmp_path):
    links = trans_links(trans_cooler, top_n=10)
    bins = trans_cooler.bins()[:][["chrom", "start", "end"]]
    track = bins.assign(E1=np.where(bins["start"] // 5_000_000 % 2 == 0, 1.0, -1.0))
    out = tmp_path / "circos.png"
    ax = plot_circos(CHROMSIZES, links, track=track, out=str(out), title="test")
    assert out.stat().st_size > 0
    texts = {t.get_text() for t in ax.texts}
    assert {"1", "2", "3"} <= texts
    # 3 cromosomi + 45 barre della traccia + 10 corde
    assert len(ax.patches) == 3 + len(track) + 10


def test_plot_circos_highlight_draws_matching_links_on_top(trans_cooler):
    links = trans_links(trans_cooler, top_n=10)
    ax = plot_circos(CHROMSIZES, links, highlight="chr3:7Mb-8Mb")
    chords = [p for p in ax.patches if isinstance(p, matplotlib.patches.PathPatch)]
    green = [p for p in chords if p.get_zorder() == 3.5]
    n_hit = ((links.chrom1 == "chr3") & (links.start1 == 7_000_000)).sum() + \
            ((links.chrom2 == "chr3") & (links.start2 == 7_000_000)).sum()
    assert len(chords) == 10 and len(green) == n_hit >= 1


def test_cis_links_finds_long_range_hotspot(trans_cooler):
    links = cis_links(trans_cooler, min_distance=10_000_000, top_n=5)
    assert (links["chrom1"] == links["chrom2"]).all()
    assert ((links["start2"] - links["start1"]) >= 10_000_000).all()
    top = links.iloc[0]
    chrom, x, y = CIS_HOTSPOT
    assert (top.chrom1, top.start1, top.start2) == (chrom, x * BINSIZE, y * BINSIZE)
    # dopo il filtro delle estremità, la diagonale a 13 Mb di un cromosoma da
    # 20 bin ha solo 3 celle e l'hotspot stesso alza la media: resta comunque primo
    assert top.oe > 1.5 * links.iloc[1].oe


def test_cis_links_respects_min_distance(trans_cooler):
    assert cis_links(trans_cooler, min_distance=25_000_000).empty  # nessun cromosoma è così lungo


def test_plot_circos_colors_cis_and_trans(trans_cooler):
    links = pd.concat([trans_links(trans_cooler, top_n=3), cis_links(trans_cooler, top_n=2)])
    ax = plot_circos(CHROMSIZES, links, link_color="#ff00ff", cis_color="#ff8000")
    chords = [p for p in ax.patches if isinstance(p, matplotlib.patches.PathPatch)]
    colors = [matplotlib.colors.to_hex(p.get_edgecolor()) for p in chords]
    assert colors.count("#ff00ff") == 3 and colors.count("#ff8000") == 2


def test_filters_drop_edges_and_sparse_chromosomes(trans_cooler):
    from clew.circos import _usable_bins

    usable = _usable_bins(trans_cooler, None, edge_bins=2, min_valid_frac=0.5, max_bins=1000)
    # le estremità (2 bin per lato) non sono utilizzabili, il resto sì
    assert not usable["chr1"][:2].any() and not usable["chr1"][-2:].any()
    assert usable["chr1"][2:-2].all()

    links = trans_links(trans_cooler, top_n=1000)
    for c, s0 in ((links.chrom1, links.start1), (links.chrom2, links.start2)):
        size = CHROMSIZES[c].to_numpy()
        assert ((s0 >= 2 * BINSIZE) & (s0 < size - 2 * BINSIZE)).all()


def test_chromosome_with_few_valid_bins_is_excluded(tmp_path):
    """Un cromosoma quasi tutto filtrato (come chrY in una linea femminile) non dà corde."""
    sizes = pd.Series({"chr1": 20_000_000, "chrY": 20_000_000})
    bins = cooler.binnify(sizes, BINSIZE)
    i, j = np.triu_indices(len(bins))
    counts = np.where(bins["chrom"].to_numpy()[i] == bins["chrom"].to_numpy()[j], 1000 // (j - i + 1), 20)
    y = np.flatnonzero(bins["chrom"].to_numpy() == "chrY")
    dead = np.isin(i, y[:15]) | np.isin(j, y[:15])   # 15 bin su 20 di chrY senza letture
    counts[dead] = 0
    counts[(i == 10) & (j == y[17])] = 500            # e un falso hotspot su un bin rimasto
    path = str(tmp_path / "y.cool")
    cooler.create_cooler(path, bins, pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": counts}),
                         ordered=True)
    clr = cooler.Cooler(path)
    cooler.balance_cooler(clr, mad_max=0, store=True)

    assert trans_links(clr).empty
    assert not trans_links(clr, min_valid_frac=0).empty  # senza filtro l'hotspot compare
