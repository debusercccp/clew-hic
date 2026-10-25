"""Registra file .mcool in un server HiGlass locale e crea una vista da aprire nel browser.

I file non vengono copiati: ``data/`` è montata nella cartella media del server
(``/data/media/hic``) e ogni .mcool viene registrato per percorso con l'API
``link_tile``. Solo le dimensioni dei cromosomi (pochi kB) vengono caricate via
upload, per le etichette degli assi. Rilanciare lo script non crea duplicati:
gli identificativi derivano dal nome del file e dal contenuto.

Uso (dopo ``make higlass``)::

    docker compose run --rm worker python -m clew.ingest_higlass \\
        data/4DNFIXP4QG5B.mcool --region chr2:10Mb-20Mb

Con più file crea viste affiancate con zoom e posizione sincronizzati.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

import cooler
import requests

from .hic_io import list_resolutions

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
MEDIA_PREFIX = "hic"  # dove docker-compose monta data/ dentro /data/media
FALL_COLORS = ["white", "rgba(245,166,35,1.0)", "rgba(208,2,27,1.0)", "black"]


def tileset_uid(media_rel: str, prefix: str = MEDIA_PREFIX) -> str:
    """Identificativo stabile del tileset dal percorso dentro data/.

    ``hic/4DNFIXP4QG5B.mcool`` -> ``4DNFIXP4QG5B``; ``hic/demo/x.mcool`` -> ``demo_x``:
    file con lo stesso nome in cartelle diverse non collidono.
    """
    rel = media_rel.removeprefix(prefix + "/")
    rel = rel[: -len(".mcool")] if rel.endswith(".mcool") else rel
    return re.sub(r"[^A-Za-z0-9_-]", "_", rel)


def media_path(path: Path, data_dir: Path = DATA_DIR, prefix: str = MEDIA_PREFIX) -> str:
    """Percorso del file visto dal server HiGlass, relativo alla sua cartella media.

    Raises:
        ValueError: se il file non è dentro ``data/`` (il server non lo vedrebbe).
    """
    try:
        rel = path.resolve().relative_to(data_dir.resolve())
    except ValueError:
        raise ValueError(f"{path} non è dentro {data_dir}: HiGlass vede solo data/") from None
    return f"{prefix}/{rel.as_posix()}"


def chromsizes_tsv(clr: cooler.Cooler) -> str:
    """Dimensioni dei cromosomi nell'ordine del file (lo stesso usato dalla heatmap)."""
    return "".join(f"{name}\t{size}\n" for name, size in clr.chromsizes.items())


def genome_domain(chromsizes: Sequence[tuple[str, int]], region: Optional[str]) -> list[int]:
    """Converte una regione in coordinate assolute HiGlass (cromosomi concatenati)."""
    sizes = dict(chromsizes)
    total = sum(sizes.values())
    if not region:
        return [0, total]
    chrom, start, end = cooler.util.parse_region(region, sizes)
    offset = 0
    for name, size in chromsizes:
        if name == chrom:
            break
        offset += size
    return [offset + start, offset + end]


def build_viewconf(
    tilesets: Sequence[tuple[str, str]],
    chromsizes_uid: str,
    domain: Sequence[int],
    server: str = "/api/v1",
) -> dict:
    """Viewconf con una heatmap per tileset, affiancate e sincronizzate.

    Args:
        tilesets: coppie (uid, nome mostrato).
        chromsizes_uid: tileset delle dimensioni dei cromosomi (etichette degli assi).
        domain: [inizio, fine] in coordinate assolute.
        server: URL dell'API vista dal browser.
    """
    width = 12 // len(tilesets)
    center = (domain[0] + domain[1]) / 2
    views = []
    for i, (uid, name) in enumerate(tilesets):
        vid = f"view{i}"
        labels = dict(server=server, tilesetUid=chromsizes_uid, options={"color": "#808080"})
        views.append({
            "uid": vid,
            "initialXDomain": list(domain),
            "initialYDomain": list(domain),
            "layout": {"i": vid, "x": i * width, "y": 0, "w": width, "h": 12,
                       "moved": False, "static": False},
            "tracks": {
                "top": [{"type": "horizontal-chromosome-labels", "uid": f"{vid}-top",
                         "height": 30, **labels}],
                "left": [{"type": "vertical-chromosome-labels", "uid": f"{vid}-left",
                          "width": 30, **labels}],
                "center": [{
                    "type": "combined", "uid": f"{vid}-center",
                    "contents": [{
                        "type": "heatmap", "uid": f"{vid}-heatmap",
                        "server": server, "tilesetUid": uid,
                        "options": {"name": name, "colorRange": FALL_COLORS,
                                    "labelPosition": "bottomRight",
                                    "colorbarPosition": "topRight",
                                    "showMousePosition": True},
                    }],
                }],
                "right": [], "bottom": [], "whole": [], "gallery": [],
            },
        })

    vids = [v["uid"] for v in views]

    def lock(name: str) -> dict:
        if len(vids) < 2:
            return {"locksByViewUid": {}, "locksDict": {}}
        return {"locksByViewUid": {v: name for v in vids},
                "locksDict": {name: {"uid": name, **{v: [center, center, 1] for v in vids}}}}

    return {
        "editable": True,
        "zoomFixed": False,
        "trackSourceServers": [server],
        "exportViewUrl": f"{server}/viewconfs/",
        "views": views,
        "zoomLocks": lock("zoom"),
        "locationLocks": lock("location"),
        "valueScaleLocks": {"locksByViewUid": {}, "locksDict": {}},
    }


class HiGlassClient:
    """Client minimo per l'API REST di higlass-server."""

    def __init__(self, url: str, user: str, password: str, timeout: float = 30):
        self.api = url.rstrip("/") + "/api/v1"
        self.session = requests.Session()
        self.session.auth = (user, password)
        self.timeout = timeout

    def exists(self, uid: str) -> bool:
        r = self.session.get(f"{self.api}/tileset_info/", params={"d": uid}, timeout=self.timeout)
        r.raise_for_status()
        return "error" not in r.json().get(uid, {"error": True})

    def link(self, filepath: str, uid: str, filetype: str, datatype: str, assembly: str) -> None:
        """Registra un file già presente nella cartella media del server."""
        r = self.session.post(f"{self.api}/link_tile/", timeout=self.timeout, json={
            "filepath": filepath, "uuid": uid, "filetype": filetype,
            "datatype": datatype, "coordSystem": assembly,
        })
        if r.status_code != 201:
            raise RuntimeError(f"link_tile {filepath}: {r.status_code} {r.text}")

    def upload(self, content: str, filename: str, uid: str, filetype: str,
               datatype: str, assembly: str) -> None:
        """Carica un file piccolo (es. chromsizes) nel server."""
        r = self.session.post(
            f"{self.api}/tilesets/", timeout=self.timeout,
            files={"datafile": (filename, content.encode())},
            data={"uid": uid, "filetype": filetype, "datatype": datatype,
                  "coordSystem": assembly, "name": filename},
        )
        if r.status_code != 201:
            raise RuntimeError(f"upload {filename}: {r.status_code} {r.text}")

    def save_viewconf(self, viewconf: dict, uid: str) -> None:
        r = self.session.post(f"{self.api}/viewconfs/", timeout=self.timeout,
                              json={"uid": uid, "viewconf": viewconf})
        if r.status_code != 200 and "already exists" not in r.text:
            raise RuntimeError(f"viewconf: {r.status_code} {r.text}")


def ingest(
    paths: Sequence[Path],
    client: HiGlassClient,
    assembly: str,
    region: Optional[str] = None,
    public_url: str = "http://127.0.0.1:8989",
) -> str:
    """Registra i file (se serve), crea la vista e ritorna l'URL da aprire."""
    tilesets, chromsizes = [], None
    for path in paths:
        resolutions = list_resolutions(path, verbose=False)
        clr = cooler.Cooler(f"{path}::/resolutions/{resolutions[0]}")
        if "weight" not in clr.bins().columns:
            print(f"attenzione: {path.name} non è bilanciato a {resolutions[0]} bp")

        tsv = chromsizes_tsv(clr)
        if chromsizes is not None and tsv != chromsizes:
            raise ValueError(f"{path.name} ha cromosomi diversi dagli altri file: non confrontabili")
        chromsizes = tsv

        filepath = media_path(path)  # prima di tutto: il server deve poter vedere il file
        uid = tileset_uid(filepath)
        if client.exists(uid):
            print(f"{path.name}: già registrato ({uid})")
        else:
            client.link(filepath, uid, "cooler", "matrix", assembly)
            print(f"{path.name}: registrato ({uid}), risoluzioni {resolutions[0]}-{resolutions[-1]} bp")
        tilesets.append((uid, path.name))

    cs_uid = f"chromsizes-{assembly}-{hashlib.sha1(chromsizes.encode()).hexdigest()[:8]}"
    if not client.exists(cs_uid):
        client.upload(chromsizes, f"{assembly}.chrom.sizes", cs_uid,
                      "chromsizes-tsv", "chromsizes", assembly)

    pairs = [tuple(line.split("\t")) for line in chromsizes.splitlines()]
    domain = genome_domain([(c, int(s)) for c, s in pairs], region)
    viewconf = build_viewconf(tilesets, cs_uid, domain)
    vc_uid = "view-" + hashlib.sha1(json.dumps(viewconf, sort_keys=True).encode()).hexdigest()[:12]
    client.save_viewconf(viewconf, vc_uid)
    return f"{public_url.rstrip('/')}/app/?config={vc_uid}"


def _main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(prog="python -m clew.ingest_higlass",
                                     description=__doc__.split("\n")[0])
    parser.add_argument("paths", nargs="+", type=Path, help="file .mcool dentro data/")
    parser.add_argument("--assembly", default="hg38", help="coordSystem (hg38, mm10...)")
    parser.add_argument("--region", help="vista iniziale, es. chr2:10Mb-20Mb (default: genoma)")
    parser.add_argument("--url", default=os.environ.get("HIGLASS_URL", "http://127.0.0.1:8989"),
                        help="API del server (dal container worker: http://higlass)")
    parser.add_argument("--public-url",
                        default=os.environ.get("HIGLASS_PUBLIC_URL", "http://127.0.0.1:8989"),
                        help="indirizzo del server visto dal browser")
    args = parser.parse_args(argv)

    user, password = os.environ.get("HIGLASS_USER", "clew"), os.environ.get("HIGLASS_PASSWORD")
    if not password:
        parser.error("HIGLASS_PASSWORD non impostata (la crea `make env` in .env)")

    try:
        url = ingest(args.paths, HiGlassClient(args.url, user, password),
                     args.assembly, args.region, args.public_url)
    except (ValueError, RuntimeError, requests.ConnectionError) as e:
        parser.exit(1, f"errore: {e}\n")
    print(f"\napri nel browser:\n  {url}")


if __name__ == "__main__":
    _main()
