"""
A tiny synthetic st-qc output with planted spatial patterns, in two layouts:

  <out>/linked/spatial/<sample_id>.h5ad, <out>/linked/metadata/samples.csv
  <out>/flat/...   the same files side by side (no links)

Samples (each follows the AnnData contract after st-qc: X log-normalized
float32 CSR, layers["counts"], var["highly_variable"], uns["pennsieve"]):

  vis_a, vis_b   Visium-like hex grid (24 x 20 = 480 spots, full-resolution pixel
                 coordinates, a lowres "H&E" image and scale factors)
  cells_a        Xenium-like cells (600 random positions in microns, no image)

Planted genes (counts ~ Poisson with a spatially varying rate):
  GRAD       a left-to-right gradient
  BLOB       a disc of high expression
  STRIPE     a band across the tissue
  PATCHY     two discs (vis_a and cells_a only; flat in vis_b)
Noise genes N000..N199 and CONST are Poisson with a constant rate (no pattern;
they make up most of each library, as unpatterned genes do in real data, so
library-size normalization doesn't induce a pattern in them); RARE is detected
in only 2 observations (not tested). N199 is not an HVG.

Usage: python -m tests.make_fixture <out> [--layout linked|flat|both] [--raw]
  --raw writes st-prep-like files (raw counts in X, no 'normalized' flag).
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

PATTERNED = ["GRAD", "BLOB", "STRIPE", "PATCHY"]
NOISE = [f"N{i:03d}" for i in range(200)]
GENES = PATTERNED + NOISE + ["CONST", "RARE"]


def hex_grid(rows: int, cols: int, spacing: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rr, cc = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    rr, cc = rr.ravel(), cc.ravel()
    x = cc * spacing + (rr % 2) * spacing / 2 + 300.0
    y = rr * spacing * np.sqrt(3) / 2 + 300.0
    return np.column_stack([x, y]), rr, cc * 2 + (rr % 2)


def rates(u: np.ndarray, v: np.ndarray, sample: str) -> dict[str, np.ndarray]:
    """u, v in [0, 1]: relative position in the tissue."""
    r = {}
    r["GRAD"] = 0.5 + 8.0 * u
    if sample == "vis_b":
        r["BLOB"] = 0.5 + 10.0 * (((u - 0.3) ** 2 + (v - 0.6) ** 2) < 0.04)
    else:
        r["BLOB"] = 0.5 + 10.0 * (((u - 0.6) ** 2 + (v - 0.4) ** 2) < 0.04)
    r["STRIPE"] = 0.5 + 8.0 * (np.abs(v - 0.5) < 0.12)
    if sample == "vis_b":
        r["PATCHY"] = np.full(u.shape, 3.0)
    else:
        r["PATCHY"] = 0.5 + 9.0 * ((((u - 0.2) ** 2 + (v - 0.2) ** 2) < 0.02) |
                                   (((u - 0.8) ** 2 + (v - 0.8) ** 2) < 0.02))
    for g in NOISE:
        r[g] = np.full(u.shape, 2.0)
    r["CONST"] = np.full(u.shape, 5.0)
    r["RARE"] = np.zeros(u.shape)
    return r


def make_sample(sample_id: str, rng: np.random.Generator, raw: bool):
    import anndata as ad
    import scanpy as sc

    if sample_id.startswith("vis"):
        coords, arow, acol = hex_grid(24, 20, 150.0)
        kind, platform, unit, grid, diameter = "spot", "visium", "pixel", "hex", 100.0
    else:
        coords = rng.uniform(0, 400, size=(600, 2))
        arow = acol = None
        kind, platform, unit, grid, diameter = "cell", "xenium", "micron", "none", 10.0
    n = coords.shape[0]
    span = coords.max(axis=0) - coords.min(axis=0)
    uv = (coords - coords.min(axis=0)) / span
    r = rates(uv[:, 0], uv[:, 1], sample_id)
    counts = np.column_stack([rng.poisson(r[g]) for g in GENES]).astype(np.float32)
    rare = np.flatnonzero(rng.permutation(n) < 2)
    counts[rare, GENES.index("RARE")] = 3
    counts = sp.csr_matrix(counts)
    obs = pd.DataFrame(index=[f"{sample_id}:obs{i:04d}" for i in range(n)])
    obs["sample_id"] = pd.Categorical([sample_id] * n)
    obs["x"], obs["y"] = coords[:, 0], coords[:, 1]
    obs["in_tissue"] = 1
    if arow is not None:
        obs["array_row"], obs["array_col"] = arow, acol
    var = pd.DataFrame(index=GENES)
    var["gene_ids"] = [f"ENSG{i:011d}" for i in range(len(GENES))]
    var["feature_types"] = "Gene Expression"
    var["mt"] = False
    a = ad.AnnData(X=counts, obs=obs, var=var)
    a.obsm["spatial"] = coords.astype(np.float64)
    pen = {"platform": platform, "kind": kind, "coordinate_unit": unit, "spot_diameter": diameter, "grid": grid,
           "source": f"{sample_id}/", "prep_version": "0.1.0"}
    if platform == "visium":
        scalef = 0.1
        h = int(coords[:, 1].max() * scalef) + 40
        w = int(coords[:, 0].max() * scalef) + 40
        yy, xx = np.mgrid[0:h, 0:w]
        img = np.stack([200 + 30 * np.sin(xx / 7.0), 150 + 40 * np.cos(yy / 9.0), 200 + 0 * xx], axis=-1)
        a.uns["spatial"] = {sample_id: {
            "images": {"lowres": np.clip(img, 0, 255).astype(np.uint8)},
            "scalefactors": {"tissue_lowres_scalef": scalef, "tissue_hires_scalef": 0.2,
                             "spot_diameter_fullres": diameter, "fiducial_diameter_fullres": diameter * 1.6}}}
    if not raw:
        a.layers["counts"] = a.X.copy()
        target = float(np.median(np.asarray(a.X.sum(axis=1)).ravel()))
        sc.pp.normalize_total(a, target_sum=target)
        sc.pp.log1p(a)
        a.X = sp.csr_matrix(a.X, dtype=np.float32)
        hv = np.ones(len(GENES), dtype=bool)
        hv[GENES.index("N199")] = False          # one gene not an HVG, so GENES=hvg and all differ
        a.var["highly_variable"] = hv
        pen.update(normalized=True, target_sum=target)
    a.uns["pennsieve"] = pen
    return a


def write(out: Path, layout: str = "both", raw: bool = False, seed: int = 7) -> None:
    rng = np.random.default_rng(seed)
    samples = {sid: make_sample(sid, rng, raw) for sid in ("vis_a", "vis_b", "cells_a")}
    table = pd.DataFrame({
        "sample_id": list(samples), "path": [f"{s}/" for s in samples],
        "platform": [a.uns["pennsieve"]["platform"] for a in samples.values()],
        "kind": [a.uns["pennsieve"]["kind"] for a in samples.values()],
        "n_obs": [a.n_obs for a in samples.values()], "n_genes": [a.n_vars for a in samples.values()],
        "n_obs_in_tissue": [a.n_obs for a in samples.values()],
        "condition": ["A", "B", "A"]})
    layouts = ["linked", "flat"] if layout == "both" else [layout]
    for lay in layouts:
        root = out / lay
        if root.exists():
            shutil.rmtree(root)
        sp_dir = root / "spatial" if lay == "linked" else root
        md_dir = root / "metadata" if lay == "linked" else root
        sp_dir.mkdir(parents=True, exist_ok=True)
        md_dir.mkdir(parents=True, exist_ok=True)
        for sid, a in samples.items():
            a.write_h5ad(sp_dir / f"{sid}.h5ad")
        table.to_csv(md_dir / "samples.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    ap.add_argument("--raw", action="store_true")
    args = ap.parse_args()
    write(args.out, args.layout, args.raw)


if __name__ == "__main__":
    main()
