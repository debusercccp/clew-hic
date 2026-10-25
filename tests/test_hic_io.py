import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from clew import hic_io
from tests.conftest import RESOLUTIONS


def test_list_resolutions_mcool(synthetic_mcool, capsys):
    assert hic_io.list_resolutions(synthetic_mcool) == RESOLUTIONS
    assert "10 kb, 20 kb, 40 kb" in capsys.readouterr().out


def test_list_resolutions_cool(synthetic_cool):
    assert hic_io.list_resolutions(synthetic_cool, verbose=False) == [10_000]


def test_hic_is_rejected_with_hint(tmp_path):
    with pytest.raises(ValueError, match="hictk convert"):
        hic_io.list_resolutions(tmp_path / "x.hic")


def test_open_cooler_unknown_resolution(synthetic_mcool):
    with pytest.raises(ValueError, match="disponibili"):
        hic_io.open_cooler(synthetic_mcool, 5_000)


def test_fetch_region_shape_and_symmetry(synthetic_mcool):
    mat = hic_io.fetch_region(synthetic_mcool, 10_000, "chr1:100,000-500,000")
    assert mat.shape == (40, 40)
    assert np.allclose(mat, mat.T, equal_nan=True)
    assert np.isfinite(mat[5:-5, 5:-5]).all()


def test_fetch_region_raw_counts(synthetic_mcool):
    raw = hic_io.fetch_region(synthetic_mcool, 10_000, "chr1:0-100000", balance=False)
    assert raw[0, 0] == round(2000 * 4 * 1.5)


def test_fetch_region_max_bins_guard(synthetic_mcool):
    with pytest.raises(ValueError, match="limite"):
        hic_io.fetch_region(synthetic_mcool, 10_000, "chr1", max_bins=100)


def test_fetch_region_requires_weights(synthetic_cool):
    with pytest.raises(ValueError, match="non è bilanciato"):
        hic_io.fetch_region(synthetic_cool, 10_000, "chr1:0-100000")


def test_plot_matrix_handles_zeros_and_nans(tmp_path):
    mat = np.random.default_rng(0).random((20, 20))
    mat[0, :] = np.nan
    mat[1, 1] = 0
    out = tmp_path / "m.png"
    ax = hic_io.plot_matrix(mat, "chr2:10,000,000-10,200,000", out=out)
    assert out.stat().st_size > 0
    assert ax.get_xlim() == pytest.approx((10.0, 10.2))
    assert ax.get_xlabel() == "chr2 (Mb)"


def test_plot_pair_shares_color_scale(synthetic_mcool, tmp_path):
    fig = hic_io.plot_pair(synthetic_mcool, 20_000, "chr1:0-1,000,000", "chr2:0-1,000,000",
                           out=tmp_path / "pair.png")
    im_a, im_b = fig.axes[0].images[0], fig.axes[1].images[0]
    assert im_a.norm.vmin == im_b.norm.vmin and im_a.norm.vmax == im_b.norm.vmax
