"""The synthetic fixture through the processor (no Docker), both input layouts."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import main as M
from processor import svg as S
from processor.graph import build_graph
from tests import make_fixture

PLANTED = {"vis_a": ["GRAD", "BLOB", "STRIPE", "PATCHY"], "vis_b": ["GRAD", "BLOB", "STRIPE"],
           "cells_a": ["GRAD", "BLOB", "STRIPE", "PATCHY"]}


@pytest.fixture(scope="session")
def fixture_dir(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("fixture")
    make_fixture.write(out, "both")
    make_fixture.write(out / "raw", "linked", raw=True)
    return out


def run(monkeypatch, input_dir: Path, output_dir: Path, **params) -> int:
    for k in M.PARAMETERS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("THREADS", "1")
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    return M.main()


@pytest.fixture(scope="session")
def default_run(fixture_dir, tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("out-default")
    mp = pytest.MonkeyPatch()
    try:
        assert run(mp, fixture_dir / "linked", out) == 0
    finally:
        mp.undo()
    return out


def test_outputs_and_columns(default_run):
    names = {p.name for p in default_run.iterdir()}
    for f in ("svg.csv", "svg_consensus.csv", "summary.json", "provenance.json", "morans_distribution.png",
              "consensus_top.png", "vis_a_svg_top.png", "vis_b_svg_top.png", "cells_a_svg_top.png"):
        assert f in names, f
    svg = pd.read_csv(default_run / "svg.csv")
    assert list(svg.columns) == ["sample_id", "gene", "morans_i", "pval", "pval_adj", "significant", "rank"]
    cons = pd.read_csv(default_run / "svg_consensus.csv")
    assert list(cons.columns) == ["gene", "n_samples", "n_significant", "mean_morans_i", "median_rank", "rank"]
    prov = json.loads((default_run / "provenance.json").read_text())
    assert prov["app"] == "st-svg" and prov["version"] == "0.1.0" and prov["n_ok"] == 3
    assert prov["tools"]["squidpy"]
    summary = json.loads((default_run / "summary.json").read_text())
    assert summary["n_samples"] == 3 and len(summary["consensus_top"]) == 12


def test_planted_genes_found(default_run):
    svg = pd.read_csv(default_run / "svg.csv")
    for sid, genes in PLANTED.items():
        d = svg[svg["sample_id"] == sid].set_index("gene")
        for g in genes:
            assert bool(d.loc[g, "significant"]), (sid, g)
            assert d.loc[g, "rank"] <= len(genes), (sid, g)
            assert d.loc[g, "morans_i"] > 0.3, (sid, g)
        # noise: few false positives
        noise = d[d.index.str.match(r"^N\d{3}$")]
        assert noise["significant"].mean() <= 0.03, sid
        assert d.index.is_unique
        assert sorted(d["rank"]) == list(range(1, len(d) + 1))
    vb = svg[svg["sample_id"] == "vis_b"].set_index("gene")
    assert not bool(vb.loc["PATCHY", "significant"])        # flat in vis_b
    assert "RARE" not in set(svg["gene"])                    # detected in 2 observations: not tested
    assert "N199" not in set(svg["gene"])                    # not an HVG
    cons = pd.read_csv(default_run / "svg_consensus.csv").set_index("gene")
    assert set(cons.index[:3]) == {"GRAD", "BLOB", "STRIPE"}
    assert cons.loc["PATCHY", "n_significant"] == 2 and cons.loc["PATCHY", "rank"] == 4
    assert (cons.loc[["GRAD", "BLOB", "STRIPE"], "n_significant"] == 3).all()


def test_flat_layout_same_result(default_run, fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir / "flat", tmp_path) == 0
    a = pd.read_csv(default_run / "svg.csv")
    b = pd.read_csv(tmp_path / "svg.csv")
    pd.testing.assert_frame_equal(a, b)


def test_mode_both_and_all_genes(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir / "linked", tmp_path, MODE="both", GENES="all", TOP_N=4) == 0
    svg = pd.read_csv(tmp_path / "svg.csv")
    assert {"gearys_c", "gearys_pval", "gearys_pval_adj"} <= set(svg.columns)
    assert "N199" in set(svg["gene"])
    d = svg[svg["sample_id"] == "vis_a"].set_index("gene")
    for g in PLANTED["vis_a"]:
        assert d.loc[g, "gearys_c"] < 0.7 and d.loc[g, "gearys_pval_adj"] < 0.05


def test_gene_list_and_permutations(fixture_dir, tmp_path, monkeypatch):
    rc = run(monkeypatch, fixture_dir / "linked", tmp_path, GENES="GRAD, BLOB,N000,NOPE", N_PERMS=99)
    assert rc == 0
    svg = pd.read_csv(tmp_path / "svg.csv")
    assert set(svg["gene"]) == {"GRAD", "BLOB", "N000"}
    planted = svg[svg["gene"].isin(["GRAD", "BLOB"])]
    assert np.allclose(planted["pval"], 0.01)                 # (0 + 1) / (99 + 1)
    assert planted["significant"].all()
    assert (svg["pval"] > 0).all() and (svg["pval"] <= 1).all()


def test_raw_counts_normalized_here(fixture_dir, tmp_path, monkeypatch):
    assert run(monkeypatch, fixture_dir / "raw" / "linked", tmp_path) == 0
    svg = pd.read_csv(tmp_path / "svg.csv")
    assert "N199" in set(svg["gene"])                         # no var['highly_variable']: <= 2000 genes, all used
    d = svg[svg["sample_id"] == "vis_a"].set_index("gene")
    assert all(bool(d.loc[g, "significant"]) for g in PLANTED["vis_a"])
    prov = json.loads((tmp_path / "provenance.json").read_text())
    assert prov["samples"][0]["normalization"].startswith("raw counts")


def test_bad_sample_skipped_and_failures(fixture_dir, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    (inp / "spatial").mkdir(parents=True)
    os.symlink(fixture_dir / "linked" / "spatial" / "vis_a.h5ad", inp / "spatial" / "vis_a.h5ad")
    (inp / "spatial" / "broken.h5ad").write_bytes(b"not hdf5")
    assert run(monkeypatch, inp, tmp_path / "out") == 0
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    assert prov["n_ok"] == 1 and prov["n_failed"] == 1
    (inp / "spatial" / "vis_a.h5ad").unlink()
    assert run(monkeypatch, inp, tmp_path / "out2") == 1     # every sample failed
    empty = tmp_path / "empty"
    empty.mkdir()
    assert run(monkeypatch, empty, tmp_path / "out3") == 1   # no input
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "o4", MODE="geary") == 2
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "o5", NEIGHBOURS="many") == 2
    assert run(monkeypatch, fixture_dir / "linked", tmp_path / "o6", FDR="1.5") == 2


def test_knn_graph_and_cap():
    rng = np.random.default_rng(0)
    pts = np.vstack([rng.uniform(0, 100, (200, 2)), rng.uniform(1000, 1100, (50, 2))])
    pen = {"grid": "none", "spot_diameter": 10.0}
    g, info = build_graph(pts, pen, 1, 6, 3.0)
    assert info["mode"] == "knn" and (g != g.T).nnz == 0 and g.diagonal().sum() == 0
    # no edge between the two clusters 900+ apart (cap 30)
    assert g[:200, 200:].nnz == 0
    g2, info2 = build_graph(pts, pen, 1, "delaunay", 3.0)
    assert info2["mode"] == "delaunay" and g2[:200, 200:].nnz == 0
    g3, _ = build_graph(pts, pen, 1, "delaunay", 0.0)         # no cap: Delaunay joins the clusters
    assert g3[:200, 200:].nnz > 0


def test_bh_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests

    p = np.random.default_rng(1).uniform(size=500) ** 3
    assert np.allclose(S.bh_adjust(p), multipletests(p, method="fdr_bh")[1])
    assert np.isnan(S.bh_adjust(np.array([np.nan, 0.1]))[0])


def test_negative_autocorrelation_not_significant():
    """A checkerboard (negative autocorrelation) must get a large one-sided p-value."""
    rows, cols = 20, 20
    rr, cc = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    coords = np.column_stack([cc.ravel(), rr.ravel()]).astype(float)
    g, _ = build_graph(coords, {"grid": "none", "spot_diameter": 1.0}, 1, 4, 1.5)
    checker = ((rr + cc) % 2).ravel().astype(np.float32) + 0.01 * np.random.default_rng(0).normal(size=rows * cols)
    smooth = (cc.ravel() / cols).astype(np.float32)
    X = sp.csr_matrix(np.column_stack([checker, smooth]).astype(np.float32) + 1.0)
    for n_perms in (0, 49):
        df = S.test_sample("t", X, g, np.array([0, 1]), np.array(["CHECK", "SMOOTH"]), "moran", n_perms, 1, 0.05)
        d = df.set_index("gene")
        assert d.loc["CHECK", "morans_i"] < -0.5 and d.loc["CHECK", "pval"] > 0.9
        assert not d.loc["CHECK", "significant"]
        assert d.loc["SMOOTH", "significant"]


def test_lambda_handler(fixture_dir, tmp_path, monkeypatch):
    from processor import handler

    for k in M.PARAMETERS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("TOP_N", "99")          # set by a previous invocation: cleared by the handler
    res = handler.handler({"inputDir": str(fixture_dir / "linked"), "outputDir": str(tmp_path),
                           "GENES": ["GRAD", "BLOB"], "THREADS": 1}, None)
    assert res["status"] == "ok"
    assert set(pd.read_csv(tmp_path / "svg.csv")["gene"]) == {"GRAD", "BLOB"}
    assert json.loads((tmp_path / "summary.json").read_text())["parameters"]["TOP_N"] == 12
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(fixture_dir / "linked"), "outputDir": str(tmp_path), "MODE": "x"}, None)
