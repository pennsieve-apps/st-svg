"""
Spatially variable genes of one sample, and the consensus over samples.

Per sample:
  1. genes: var["highly_variable"] (GENES=hvg), every gene (all), or a list;
     genes detected in fewer than MIN_DETECTED observations or constant are not
     tested (their statistic is undefined or meaningless).
  2. squidpy.gr.spatial_autocorr(mode="moran") on the log-normalized X over the
     row-normalized spatial graph; with MODE=both also mode="geary".
  3. p-values: analytic (normality assumption) when N_PERMS=0, else from
     N_PERMS permutations. Both are made one-sided for *positive*
     autocorrelation (Moran's I above its expectation -1/(n-1), Geary's C
     below 1): squidpy reports the tail the statistic falls in.
  4. Benjamini-Hochberg over the genes tested in the sample; a gene is
     significant when pval_adj < FDR and I > E[I]; rank 1 = highest I.
Consensus: per gene over the samples it was tested in: n_samples,
n_significant, mean_morans_i, median_rank; ranked by n_significant, then
mean_morans_i.
"""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm

MIN_DETECTED = 3       # observations with a nonzero value for a gene to be tested
GENE_CHUNK = 1000      # genes per spatial_autocorr call (bounds the permutation arrays)


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN kept as NaN)."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(m)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def select_genes(X: sp.csr_matrix, var_names: np.ndarray, hvg: np.ndarray | None, spec: str | list[str]
                 ) -> tuple[np.ndarray, dict[str, Any]]:
    """Indices of the genes to test, and what was chosen and dropped."""
    info: dict[str, Any] = {}
    n_vars = len(var_names)
    if spec == "all":
        cand = np.arange(n_vars)
        info["source"] = "all genes"
    elif spec == "hvg":
        if hvg is not None and hvg.any():
            cand = np.flatnonzero(hvg)
            info["source"] = "var['highly_variable']"
        else:
            cand = hvg_fallback(X)
            info["source"] = "no var['highly_variable']: scanpy highly_variable_genes(flavor='seurat', n_top_genes=2000) here"
    else:
        pos = {g: i for i, g in enumerate(var_names)}
        cand = np.array([pos[g] for g in spec if g in pos], dtype=int)
        info["source"] = "GENES list"
        info["not_found"] = [g for g in spec if g not in pos]
    detected = np.bincount(X.indices, minlength=n_vars)
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    variance = mean_sq - mean ** 2
    keep = (detected[cand] >= MIN_DETECTED) & (variance[cand] > 1e-12)
    info["n_candidates"] = int(cand.size)
    info["n_not_tested"] = int((~keep).sum())
    return cand[keep], info


def hvg_fallback(X: sp.csr_matrix, n_top: int = 2000) -> np.ndarray:
    import anndata as ad
    import scanpy as sc

    if X.shape[1] <= n_top:
        return np.arange(X.shape[1])
    a = ad.AnnData(X)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=n_top)
    return np.flatnonzero(a.var["highly_variable"].to_numpy())


def _one_sided(stat: np.ndarray, p: np.ndarray, positive_when_above: bool, centre: float,
               n_perms: int) -> np.ndarray:
    """squidpy's permutation p-value counts the tail the statistic lies in, (k + 1) / (n + 1); for a
    statistic on the negative-autocorrelation side, turn it into the positive-autocorrelation tail."""
    wrong_tail = stat <= centre if positive_when_above else stat >= centre
    q = p.copy()
    q[wrong_tail] = np.minimum(1.0, (n_perms + 2) / (n_perms + 1) - p[wrong_tail])
    return q


def autocorr(X: sp.csr_matrix, graph: sp.csr_matrix, gene_idx: np.ndarray, var_names: np.ndarray, mode: str,
             n_perms: int, seed: int) -> pd.DataFrame:
    """Moran's I or Geary's C per gene with one-sided p-values, gene order kept."""
    import anndata as ad
    import squidpy as sq

    n = X.shape[0]
    stat_col = "I" if mode == "moran" else "C"
    frames = []
    for c, start in enumerate(range(0, gene_idx.size, GENE_CHUNK)):
        idx = gene_idx[start:start + GENE_CHUNK]
        a = ad.AnnData(X[:, idx], obsp={"spatial_connectivities": graph})
        a.var_names = var_names[idx]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df = sq.gr.spatial_autocorr(a, mode=mode, transformation=True, n_perms=n_perms or None,
                                        two_tailed=False, corr_method=None, rng=seed + c, copy=True,
                                        n_jobs=1, backend="threading", show_progress_bar=False)
        df = df.reindex(a.var_names)
        stat = df[stat_col].to_numpy(dtype=float)
        if n_perms:
            p = df["pval_sim"].to_numpy(dtype=float)
            if mode == "moran":
                p = _one_sided(stat, p, True, -1.0 / (n - 1), n_perms)
            else:
                p = _one_sided(stat, p, False, 1.0, n_perms)
        else:
            # analytic: the upper (Moran) / lower (Geary) normal tail from squidpy's variance under normality,
            # with the survival function so tiny p-values don't underflow to 0 as 1 - cdf does
            se = np.sqrt(df["var_norm"].to_numpy(dtype=float))
            p = norm.sf((stat + 1.0 / (n - 1)) / se) if mode == "moran" else norm.cdf((stat - 1.0) / se)
        frames.append(pd.DataFrame({"gene": a.var_names.to_numpy(), "stat": stat, "pval": p}))
    if not frames:
        return pd.DataFrame(columns=["gene", "stat", "pval"])
    return pd.concat(frames, ignore_index=True)


def test_sample(sample_id: str, X: sp.csr_matrix, graph: sp.csr_matrix, gene_idx: np.ndarray,
                var_names: np.ndarray, mode: str, n_perms: int, seed: int, fdr: float) -> pd.DataFrame:
    """svg.csv rows of one sample, sorted by rank."""
    m = autocorr(X, graph, gene_idx, var_names, "moran", n_perms, seed)
    expected = -1.0 / (X.shape[0] - 1)
    out = pd.DataFrame({"sample_id": sample_id, "gene": m["gene"], "morans_i": m["stat"], "pval": m["pval"]})
    out["pval_adj"] = bh_adjust(out["pval"].to_numpy())
    out["significant"] = (out["pval_adj"] < fdr) & (out["morans_i"] > expected)
    if mode == "both":
        g = autocorr(X, graph, gene_idx, var_names, "geary", n_perms, seed)
        out["gearys_c"] = g["stat"].to_numpy()
        out["gearys_pval"] = g["pval"].to_numpy()
        out["gearys_pval_adj"] = bh_adjust(g["pval"].to_numpy())
    out = out.sort_values(["morans_i", "gene"], ascending=[False, True], na_position="last", kind="mergesort")
    out["rank"] = np.arange(1, len(out) + 1)
    return out.reset_index(drop=True)


def consensus(per_sample: list[pd.DataFrame]) -> pd.DataFrame:
    cols = ["gene", "n_samples", "n_significant", "mean_morans_i", "median_rank", "rank"]
    if not per_sample:
        return pd.DataFrame(columns=cols)
    df = pd.concat(per_sample, ignore_index=True)
    g = df.groupby("gene", sort=False)
    c = pd.DataFrame({
        "n_samples": g["sample_id"].nunique(),
        "n_significant": g["significant"].sum().astype(int),
        "mean_morans_i": g["morans_i"].mean(),
        "median_rank": g["rank"].median(),
    }).reset_index()
    c = c.sort_values(["n_significant", "mean_morans_i", "gene"], ascending=[False, False, True], kind="mergesort")
    c["rank"] = np.arange(1, len(c) + 1)
    return c[cols].reset_index(drop=True)
