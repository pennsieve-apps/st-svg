"""
The spatial neighbour graph of one sample, built with squidpy as the workflow
conventions say ("Spatial neighbours"):

  grid: hex (Visium)   squidpy grid neighbours, n_neighs=6, n_rings=N_RINGS
  everything else      generic coordinates: the Delaunay triangulation or the K
                       nearest neighbours (NEIGHBOURS), with edges longer than
                       MAX_NEIGHBOUR_DISTANCE x spot_diameter dropped so sparse
                       regions don't connect across gaps.

squidpy 1.8 splits `spatial_neighbors` (deprecated) into one builder per mode;
the builders are used when present, with the same arguments, else the old call.
The result is binary and symmetric, without self loops; spatial_autocorr
row-normalizes it.
"""
from __future__ import annotations

import logging
import warnings
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _binary_symmetric(adj) -> sp.csr_matrix:
    a = sp.csr_matrix(adj, dtype=np.float32)
    a = a.maximum(a.T).tocsr()
    a.setdiag(0)
    a.eliminate_zeros()
    a.sum_duplicates()
    a.data[:] = 1.0
    return a


def spot_diameter(pennsieve: dict[str, Any]) -> float | None:
    try:
        v = float(pennsieve.get("spot_diameter"))
    except (TypeError, ValueError):
        return None
    return v if v > 0 and np.isfinite(v) else None


def build_graph(coords: np.ndarray, pennsieve: dict[str, Any], n_rings: int, neighbours: str | int,
                max_neighbour_distance: float) -> tuple[sp.csr_matrix, dict[str, Any]]:
    """Binary symmetric CSR connectivities and a description of how they were built."""
    import anndata as ad
    import scanpy as sc
    import squidpy as sq

    sc.settings.verbosity = 0
    for name in ("spatialdata", "spatialdata._logging"):   # squidpy's graph builders log through it
        logging.getLogger(name).setLevel(logging.WARNING)
    n = coords.shape[0]
    grid = str(pennsieve.get("grid") or "none").lower()
    spot = spot_diameter(pennsieve)
    cap = max_neighbour_distance * spot if (spot and max_neighbour_distance > 0) else None
    adata = ad.AnnData(obs=pd.DataFrame(index=pd.Index([str(i) for i in range(n)])),
                       obsm={"spatial": np.asarray(coords, dtype=np.float64)})
    info: dict[str, Any] = {"grid": grid, "spot_diameter": spot}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=DeprecationWarning)
        if grid == "hex":
            info.update(mode="grid", n_neighs=6, n_rings=n_rings)
            if hasattr(sq.gr, "spatial_neighbors_grid"):
                sq.gr.spatial_neighbors_grid(adata, n_neighs=6, n_rings=n_rings)
            else:
                sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=n_rings)
            adj = adata.obsp["spatial_connectivities"]
        elif neighbours == "delaunay":
            info.update(mode="delaunay", radius_cap=cap)
            radius = (0.0, cap) if cap else None
            if hasattr(sq.gr, "spatial_neighbors_delaunay"):
                sq.gr.spatial_neighbors_delaunay(adata, radius=radius)
            else:
                sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True, radius=radius)
            adj = adata.obsp["spatial_connectivities"]
        else:
            k = max(1, min(int(neighbours), n - 1))
            info.update(mode="knn", n_neighs=k, radius_cap=cap)
            if hasattr(sq.gr, "spatial_neighbors_knn"):
                sq.gr.spatial_neighbors_knn(adata, n_neighs=k)
            else:
                sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=k)
            adj = sp.csr_matrix(adata.obsp["spatial_connectivities"])
            if cap:
                dist = sp.csr_matrix(adata.obsp["spatial_distances"])
                if dist.nnz:
                    dist.data = (dist.data <= cap).astype(np.float32)
                    adj = sp.csr_matrix(adj.multiply(dist))
    graph = _binary_symmetric(adj)
    degree = np.diff(graph.indptr)
    info.update(n_edges=int(graph.nnz // 2), mean_degree=round(float(degree.mean()), 3) if n else 0.0,
                n_isolated=int((degree == 0).sum()))
    return graph, info
