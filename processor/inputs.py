"""
Finding the inputs and reading one sample.

Linked layout (placement: subfolder): INPUT_DIR/spatial/<sample_id>.h5ad and
INPUT_DIR/metadata/samples.csv (optional; gives the sample order and extra
columns such as condition). Without links the same files are found by name
anywhere under INPUT_DIR (symlinks followed).

A sample is read lean: X, obs, var, obsm["spatial"], uns["pennsieve"] and, for
Visium, the lowres image and scale factors. layers["counts"], the hires image,
PCA and the rest of the file are not loaded.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger("st-svg")


@dataclass
class Inputs:
    samples: list[tuple[str, Path]]              # (sample_id, h5ad), in samples.csv order when there is one
    samples_csv: Path | None
    rows: dict[str, dict[str, Any]]              # samples.csv row per sample_id
    files_seen: list[Path]
    skipped: dict[str, str] = field(default_factory=dict)


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinks (INPUT_DIR entries are symlinks), sorted."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def _find(root: Path, port: str, predicate) -> list[Path]:
    """Files of a port: its own folder when it exists, else anywhere under INPUT_DIR."""
    folder = root / port
    base = folder if folder.is_dir() else root
    return [p for p in walk(base) if predicate(p)]


def discover(input_dir: Path) -> Inputs:
    h5ads = _find(input_dir, "spatial", lambda p: p.suffix.lower() == ".h5ad")
    csvs = _find(input_dir, "metadata", lambda p: p.name.lower() == "samples.csv")
    if len(csvs) > 1:
        log.warning("%d files named samples.csv; using %s", len(csvs), csvs[0])
    samples_csv = csvs[0] if csvs else None

    by_id: dict[str, Path] = {}
    skipped: dict[str, str] = {}
    for p in h5ads:
        if p.stem in by_id:
            skipped[str(p)] = f"a second file for sample {p.stem}"
            continue
        by_id[p.stem] = p

    rows: dict[str, dict[str, Any]] = {}
    order: list[str] = []
    if samples_csv is not None:
        try:
            table = pd.read_csv(samples_csv, dtype=str, keep_default_na=False)
            if "sample_id" in table.columns:
                for _, r in table.iterrows():
                    sid = str(r["sample_id"]).strip()
                    rows[sid] = {k: v for k, v in r.items()}
                    if sid in by_id:
                        order.append(sid)
                    else:
                        skipped[sid] = "listed in samples.csv but no <sample_id>.h5ad found"
            else:
                log.warning("%s has no sample_id column; ignored", samples_csv)
        except (OSError, ValueError) as exc:
            log.warning("cannot read %s: %s", samples_csv, exc)
    order += sorted(s for s in by_id if s not in order)
    return Inputs([(s, by_id[s]) for s in order], samples_csv, rows, h5ads, skipped)


@dataclass
class Sample:
    sample_id: str
    X: sp.csr_matrix                 # (n_obs, n_genes) float32
    obs_names: np.ndarray
    var_names: np.ndarray
    hvg: np.ndarray | None           # var["highly_variable"] when present
    coords: np.ndarray               # (n_obs, 2) float64, obsm["spatial"]
    pennsieve: dict[str, Any]
    image: np.ndarray | None = None  # lowres H&E (uint8)
    image_scale: float | None = None # tissue_lowres_scalef
    normalized_here: bool = False


def _to_py(v):
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.ndim == 0:
        return v.item()
    return v


def read_sample(path: Path, sample_id: str) -> Sample:
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        if "X" not in f or "obsm" not in f or "spatial" not in f["obsm"]:
            raise ValueError("no X or obsm['spatial'] in the file")
        X = read_elem(f["X"])
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        coords = np.asarray(read_elem(f["obsm/spatial"]), dtype=np.float64)[:, :2]
        pennsieve: dict[str, Any] = {}
        if "uns" in f and "pennsieve" in f["uns"]:
            pennsieve = {k: _to_py(v) for k, v in read_elem(f["uns/pennsieve"]).items()}
        image, scale = None, None
        if "uns" in f and "spatial" in f["uns"]:
            libs = f["uns/spatial"]
            lib = sample_id if sample_id in libs else (next(iter(libs.keys()), None))
            if lib is not None:
                g = libs[lib]
                if "images" in g and "lowres" in g["images"] and "scalefactors" in g \
                        and "tissue_lowres_scalef" in g["scalefactors"]:
                    image = np.asarray(g["images/lowres"][()])
                    scale = float(np.asarray(g["scalefactors/tissue_lowres_scalef"][()]))
    if sp.issparse(X):
        X = sp.csr_matrix(X, dtype=np.float32)
    else:
        X = sp.csr_matrix(np.asarray(X, dtype=np.float32))
    X.sum_duplicates()
    X.eliminate_zeros()
    if "sample_id" in obs.columns and len(obs):
        sid = str(obs["sample_id"].iloc[0])
        if sid != sample_id:
            log.warning("%s: obs['sample_id'] is %s; using the file name", sample_id, sid)
    hvg = None
    if "highly_variable" in var.columns:
        hvg = var["highly_variable"].to_numpy().astype(bool)
    return Sample(sample_id, X, obs.index.to_numpy().astype(str), var.index.to_numpy().astype(str), hvg,
                  coords, pennsieve, image, scale)


def ensure_log_normalized(s: Sample) -> str:
    """X should be log-normalized (st-qc). Raw integer counts (st-prep output) are normalized here."""
    if s.pennsieve.get("normalized") in (True, "true", "True", 1):
        return "log-normalized by st-qc"
    data = s.X.data
    sample = data[: min(data.size, 200_000)]
    if sample.size and np.allclose(sample, np.round(sample)) and sample.max() > 20:
        import anndata as ad
        import scanpy as sc

        a = ad.AnnData(s.X)
        sc.pp.normalize_total(a)
        sc.pp.log1p(a)
        s.X = sp.csr_matrix(a.X, dtype=np.float32)
        s.normalized_here = True
        return "raw counts: normalize_total (median) + log1p here"
    return "assumed log-normalized (no uns['pennsieve']['normalized'], values not counts)"
