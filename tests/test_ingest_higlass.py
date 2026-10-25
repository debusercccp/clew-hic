from pathlib import Path

import pytest

from clew import hic_io
from clew.ingest_higlass import build_viewconf, chromsizes_tsv, genome_domain, media_path, tileset_uid

CHROMS = [("chr1", 2_000_000), ("chr2", 1_600_000)]


def test_tileset_uid_is_stable_and_safe():
    assert tileset_uid("hic/4DNFIXP4QG5B.mcool") == "4DNFIXP4QG5B"
    assert tileset_uid("hic/GM 12878.v2.mcool") == "GM_12878_v2"
    # stesso nome in cartelle diverse -> uid diversi
    assert tileset_uid("hic/a/x.mcool") != tileset_uid("hic/b/x.mcool")


def test_media_path_maps_data_dir(tmp_path):
    (tmp_path / "sub").mkdir()
    f = tmp_path / "sub" / "x.mcool"
    f.touch()
    assert media_path(f, data_dir=tmp_path) == "hic/sub/x.mcool"


def test_media_path_rejects_files_outside_data(tmp_path):
    with pytest.raises(ValueError, match="vede solo data/"):
        media_path(Path("/etc/passwd"), data_dir=tmp_path)


def test_chromsizes_tsv_keeps_file_order(synthetic_mcool):
    clr = hic_io.open_cooler(synthetic_mcool, 10_000)
    assert chromsizes_tsv(clr) == "chr1\t2000000\nchr2\t1600000\n"


@pytest.mark.parametrize("region, expected", [
    (None, [0, 3_600_000]),
    ("chr1:100,000-200,000", [100_000, 200_000]),
    ("chr2:100,000-200,000", [2_100_000, 2_200_000]),
    ("chr2", [2_000_000, 3_600_000]),
])
def test_genome_domain(region, expected):
    assert genome_domain(CHROMS, region) == expected


def test_viewconf_single_view_has_no_locks():
    vc = build_viewconf([("A", "a.mcool")], "cs", [0, 10])
    assert len(vc["views"]) == 1 and vc["views"][0]["layout"]["w"] == 12
    heatmap = vc["views"][0]["tracks"]["center"][0]["contents"][0]
    assert heatmap["tilesetUid"] == "A" and heatmap["type"] == "heatmap"
    assert vc["zoomLocks"]["locksDict"] == {}


def test_viewconf_pair_is_side_by_side_and_locked():
    vc = build_viewconf([("A", "a"), ("B", "b")], "cs", [0, 10])
    (v0, v1) = vc["views"]
    assert (v0["layout"]["x"], v1["layout"]["x"]) == (0, 6)
    for kind in ("zoomLocks", "locationLocks"):
        locks = vc[kind]
        lock_id = locks["locksByViewUid"]["view0"]
        assert locks["locksByViewUid"]["view1"] == lock_id
        assert set(locks["locksDict"][lock_id]) == {"uid", "view0", "view1"}
