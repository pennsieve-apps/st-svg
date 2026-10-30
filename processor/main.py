"""
st-svg: which genes have a spatial pattern?

Reads the per-sample AnnData files of st-qc (port `spatial`, X log-normalized)
and st-prep's samples.csv (port `metadata`, optional), and writes at the root
of the single output port `svg`:

  svg.csv              sample_id, gene, morans_i, pval, pval_adj, significant, rank
                       (+ gearys_c, gearys_pval, gearys_pval_adj with MODE=both)
  svg_consensus.csv    gene, n_samples, n_significant, mean_morans_i, median_rank, rank
  summary.json         headline numbers per sample and the consensus top genes
  <sample_id>_svg_top.png, morans_distribution.png, consensus_top.png
  provenance.json

Samples are processed one at a time (read, graph, test, figure, release).
See processor/svg.py for the method.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-svg")
log.setLevel(logging.INFO)   # the Lambda runtime sets the root logger to WARNING first
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-svg"
REPO_ROOT = Path(__file__).resolve().parent.parent
GENE_LIST_RE = re.compile(r"^[^,]+(,[^,]+)*$")

DEFAULTS = {
    "GENES": "hvg",
    "MODE": "moran",
    "N_PERMS": "0",
    "N_RINGS": "1",
    "NEIGHBOURS": "delaunay",
    "MAX_NEIGHBOUR_DISTANCE": "3.0",
    "FDR": "0.05",
    "TOP_N": "12",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
SVG_COLUMNS = ["sample_id", "gene", "morans_i", "pval", "pval_adj", "significant", "rank"]
GEARY_COLUMNS = ["gearys_c", "gearys_pval", "gearys_pval_adj"]


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    genes = raw["GENES"].strip()
    if genes.lower() in ("hvg", "all"):
        p["GENES"] = genes.lower()
    else:
        names = [g.strip() for g in genes.split(",") if g.strip()]
        if not names:
            raise ParameterError("GENES must be 'hvg', 'all' or a comma-separated list of gene symbols")
        p["GENES"] = list(dict.fromkeys(names))
    p["MODE"] = raw["MODE"].strip().lower()
    if p["MODE"] not in ("moran", "both"):
        raise ParameterError(f"MODE must be 'moran' or 'both', got {raw['MODE']!r}")
    nb = raw["NEIGHBOURS"].strip().lower()
    if nb == "delaunay":
        p["NEIGHBOURS"] = "delaunay"
    else:
        try:
            p["NEIGHBOURS"] = int(nb)
        except ValueError as exc:
            raise ParameterError(f"NEIGHBOURS must be 'delaunay' or an integer k, got {raw['NEIGHBOURS']!r}") from exc
        if p["NEIGHBOURS"] < 1:
            raise ParameterError("NEIGHBOURS k must be >= 1")
    try:
        for k in ("N_PERMS", "N_RINGS", "TOP_N", "SEED", "THREADS"):
            p[k] = int(raw[k])
        p["MAX_NEIGHBOUR_DISTANCE"] = float(raw["MAX_NEIGHBOUR_DISTANCE"])
        p["FDR"] = float(raw["FDR"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["N_PERMS"] < 0 or p["N_RINGS"] < 1 or p["TOP_N"] < 0 or p["THREADS"] < 0 or p["MAX_NEIGHBOUR_DISTANCE"] < 0:
        raise ParameterError("N_PERMS, TOP_N, THREADS, MAX_NEIGHBOUR_DISTANCE must be >= 0 and N_RINGS >= 1")
    if not 0 < p["FDR"] < 1:
        raise ParameterError("FDR must be in (0, 1)")
    return p


def cgroup_cpus() -> int:
    """CPUs from the cgroup quota (v2 cpu.max, v1 cfs), else the affinity mask; never the host's count alone."""
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()[:2]
        if quota != "max":
            return max(1, math.floor(int(quota) / int(period)))
    except (OSError, ValueError):
        pass
    try:
        q = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read_text())
        per = int(Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read_text())
        if q > 0 and per > 0:
            return max(1, q // per)
    except (OSError, ValueError):
        pass
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def set_threads(threads: int) -> int:
    n = threads if threads > 0 else cgroup_cpus()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = str(n)
    if "numba" in sys.modules:
        # Already imported (a warm Lambda, or tests): numba's pool size is fixed, set the threads used.
        import numba

        numba.set_num_threads(max(1, min(n, numba.config.NUMBA_NUM_THREADS)))
    else:
        os.environ["NUMBA_NUM_THREADS"] = str(n)
    return n


def _rel(path, root: Path):
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _clean_json(obj):
    if isinstance(obj, float):
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        try:
            return _clean_json(obj.item())
        except (ValueError, AttributeError):
            return str(obj)
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def _fmt(v):
    if isinstance(v, bool) or type(v).__name__ == "bool_":
        return "true" if v else "false"
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    return v


class CsvAppender:
    """svg.csv is written sample by sample, so no table of every sample is held."""

    def __init__(self, path: Path, columns: list[str]):
        self.fh = open(path, "w", newline="")
        self.w = csv.writer(self.fh)
        self.columns = columns
        self.w.writerow(columns)

    def add(self, df) -> None:
        for row in df[self.columns].itertuples(index=False):
            self.w.writerow([_fmt(v.item() if hasattr(v, "item") else v) for v in row])

    def close(self) -> None:
        self.fh.close()


def tool_versions() -> dict:
    from importlib.metadata import PackageNotFoundError, version

    out = {"python": platform.python_version()}
    for pkg in ("squidpy", "scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib", "h5py", "numba"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            pass
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)   # bytes on macOS, KiB on Linux
    except (ImportError, OSError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 1
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    threads = set_threads(params["THREADS"])

    import numpy as np
    import pandas as pd

    from processor import figures as F
    from processor import svg as S
    from processor.graph import build_graph, spot_diameter
    from processor.inputs import discover, ensure_log_normalized, read_sample

    found = discover(input_dir)
    for name, reason in found.skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    if not found.samples:
        log.error("no <sample_id>.h5ad found under %s", input_dir)
        return 1
    genes_desc = params["GENES"] if isinstance(params["GENES"], str) else f"{len(params['GENES'])} listed"
    log.info("%d samples; genes %s, mode %s, %s p-values, FDR %g, %d threads", len(found.samples), genes_desc,
             params["MODE"], f"{params['N_PERMS']} permutations" if params["N_PERMS"] else "analytic",
             params["FDR"], threads)

    output_dir.mkdir(parents=True, exist_ok=True)
    columns = SVG_COLUMNS + (GEARY_COLUMNS if params["MODE"] == "both" else [])
    table = CsvAppender(output_dir / "svg.csv", columns)
    statuses, kept, dist, figures = [], [], {}, []
    try:
        for sid, path in found.samples:
            t = time.monotonic()
            status = {"sample_id": sid, "file": _rel(path, input_dir)}
            try:
                s = read_sample(path, sid)
                norm = ensure_log_normalized(s)
                if s.X.shape[0] < 10:
                    raise ValueError(f"only {s.X.shape[0]} observations")
                graph, ginfo = build_graph(s.coords, s.pennsieve, params["N_RINGS"], params["NEIGHBOURS"],
                                           params["MAX_NEIGHBOUR_DISTANCE"])
                t_graph = time.monotonic() - t
                idx, ginfo_genes = S.select_genes(s.X, s.var_names, s.hvg, params["GENES"])
                if ginfo_genes.get("not_found"):
                    log.warning("%s: GENES not in the sample: %s", sid, ", ".join(ginfo_genes["not_found"][:20]))
                if idx.size == 0:
                    raise ValueError("no gene to test (none selected, or all undetected/constant)")
                if params["N_PERMS"]:
                    need = math.ceil(idx.size / (params["N_PERMS"] + 1) / params["FDR"])
                    if need > 1:
                        log.warning("%s: with N_PERMS=%d the smallest p-value is %.2g; with %d genes tested, BH "
                                    "rejects only when at least %d genes reach it. Raise N_PERMS or use N_PERMS=0",
                                    sid, params["N_PERMS"], 1 / (params["N_PERMS"] + 1), idx.size, need)
                df = S.test_sample(sid, s.X, graph, idx, s.var_names, params["MODE"], params["N_PERMS"],
                                   params["SEED"], params["FDR"])
                t_test = time.monotonic() - t - t_graph
                table.add(df)
                keep_cols = ["sample_id", "gene", "morans_i", "significant", "rank"]
                kept.append(df[keep_cols].copy())
                dist[sid] = (df["morans_i"].to_numpy(dtype=float), df["significant"].to_numpy(dtype=bool))
                n_sig = int(df["significant"].sum())
                top = df.head(params["TOP_N"])
                if params["TOP_N"] > 0 and len(top):
                    try:
                        pos = {g: i for i, g in enumerate(s.var_names)}
                        cols = [pos[g] for g in top["gene"]][: F.MAX_PANELS]
                        dense = s.X[:, cols].toarray()
                        labels = [f"I = {r.morans_i:.3f}, q " + (f"= {r.pval_adj:.2g}" if r.pval_adj > 0 else "< 1e-300") + ("" if r.significant else " (n.s.)")
                                  for r in top.head(F.MAX_PANELS).itertuples()]
                        figures.append(F.spatial_top(
                            sid, s.coords, [dense[:, j] for j in range(dense.shape[1])], list(top["gene"])[: F.MAX_PANELS],
                            labels, s.image, s.image_scale, spot_diameter(s.pennsieve),
                            output_dir / f"{F.safe_name(sid)}_svg_top.png"))
                    except Exception as exc:   # a figure is a convenience
                        log.exception("%s: figure failed: %s", sid, exc)
                status.update(
                    status="ok", n_obs=int(s.X.shape[0]), n_genes=int(s.X.shape[1]),
                    platform=s.pennsieve.get("platform"), kind=s.pennsieve.get("kind"), normalization=norm,
                    graph=ginfo, genes=ginfo_genes, n_genes_tested=int(len(df)), n_significant=n_sig,
                    fraction_significant=round(n_sig / len(df), 4) if len(df) else None,
                    median_morans_i=float(np.nanmedian(df["morans_i"])) if len(df) else None,
                    top_genes=[{"gene": r.gene, "morans_i": float(r.morans_i), "pval_adj": float(r.pval_adj),
                                "significant": bool(r.significant)} for r in top.itertuples()],
                    seconds={"graph": round(t_graph, 2), "test": round(t_test, 2),
                             "total": round(time.monotonic() - t, 2)},
                    condition=found.rows.get(sid, {}).get("condition"))
                log.info("%s: %d obs, graph %s (mean degree %.2f, %d isolated), %d genes tested, %d significant "
                         "(FDR %g), top %s, %.1fs", sid, s.X.shape[0], ginfo["mode"], ginfo["mean_degree"],
                         ginfo["n_isolated"], len(df), n_sig, params["FDR"],
                         ", ".join(f"{g} {i:.2f}" for g, i in zip(top["gene"][:5], top["morans_i"][:5])),
                         time.monotonic() - t)
                del s, graph, df
            except Exception as exc:   # one bad sample is skipped, not fatal
                status.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                              seconds={"total": round(time.monotonic() - t, 2)})
                log.error("%s: failed: %s", sid, exc)
            statuses.append(status)
    finally:
        table.close()
    ok = [s for s in statuses if s["status"] == "ok"]
    if not ok:
        log.error("every sample failed")
        return 1

    cons = S.consensus(kept)
    cons_out = cons.copy()
    with open(output_dir / "svg_consensus.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(list(cons_out.columns))
        for row in cons_out.itertuples(index=False):
            w.writerow([_fmt(v.item() if hasattr(v, "item") else v) for v in row])

    t_fig = time.monotonic()
    try:
        figures.append(F.morans_distribution(dist, output_dir / "morans_distribution.png"))
        top_c = cons.head(max(params["TOP_N"], 1))
        if len(top_c):
            sids = [s["sample_id"] for s in ok]
            allrows = pd.concat(kept, ignore_index=True)
            mi = allrows.pivot_table(index="gene", columns="sample_id", values="morans_i", aggfunc="first")
            sg = allrows.pivot_table(index="gene", columns="sample_id", values="significant", aggfunc="first")
            mi = mi.reindex(index=top_c["gene"], columns=sids)
            sg = sg.reindex(index=top_c["gene"], columns=sids)
            figures.append(F.consensus_top(list(top_c["gene"]), sids, mi.to_numpy(dtype=float),
                                           sg.fillna(False).to_numpy(dtype=bool), list(top_c["n_significant"]),
                                           output_dir / "consensus_top.png"))
    except Exception as exc:
        log.exception("figures failed: %s", exc)
    t_fig = time.monotonic() - t_fig

    n_ok = len(ok)
    summary = {
        "n_samples": n_ok, "n_failed": len(statuses) - n_ok,
        "statistic": "moran" if params["MODE"] == "moran" else "moran+geary",
        "pvalues": f"{params['N_PERMS']} permutations" if params["N_PERMS"] else "analytic (normality)",
        "fdr": params["FDR"],
        "n_genes_tested": int(len(cons)),
        "n_genes_significant_any": int((cons["n_significant"] > 0).sum()),
        "n_genes_significant_all": int(((cons["n_significant"] == cons["n_samples"]) & (cons["n_samples"] == n_ok)).sum()),
        "per_sample": {s["sample_id"]: {k: s.get(k) for k in ("n_obs", "n_genes_tested", "n_significant",
                                                              "fraction_significant", "median_morans_i", "condition")}
                       | {"top_genes": [g["gene"] for g in s["top_genes"]], "graph_mean_degree": s["graph"]["mean_degree"]}
                       for s in ok},
        "consensus_top": [{"gene": r.gene, "n_samples": int(r.n_samples), "n_significant": int(r.n_significant),
                           "mean_morans_i": float(r.mean_morans_i), "median_rank": float(r.median_rank)}
                          for r in cons.head(params["TOP_N"]).itertuples()],
        "parameters": params,
        "figures": [p.name for p in figures],
    }
    write_json(output_dir / "summary.json", summary)
    write_json(output_dir / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "threads": threads,
        "method": {
            "graph": "squidpy spatial neighbours per conventions.md: hex grid -> grid n_neighs=6, n_rings=N_RINGS; "
                     "else generic Delaunay or kNN capped at MAX_NEIGHBOUR_DISTANCE x spot_diameter; binary, symmetric",
            "statistic": "squidpy.gr.spatial_autocorr(mode='moran'[, 'geary'], transformation=True) on log-normalized X",
            "pvalues": "analytic under normality" if not params["N_PERMS"] else f"{params['N_PERMS']} permutations (pval_sim)",
            "tail": "one-sided, positive autocorrelation (I > -1/(n-1); C < 1)",
            "correction": "Benjamini-Hochberg per sample over the genes tested",
            "significant": "pval_adj < FDR and I > -1/(n-1)",
            "min_detected": S.MIN_DETECTED,
            "consensus_rank": "n_significant desc, then mean_morans_i desc",
        },
        "inputs": {"spatial": [_rel(p, input_dir) for p in found.files_seen],
                   "samples_csv": _rel(found.samples_csv, input_dir), "skipped": found.skipped},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"figures_cohort": round(t_fig, 2), "total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(),
        "samples": statuses, "n_ok": n_ok, "n_failed": len(statuses) - n_ok,
        "outputs": sorted(p.name for p in output_dir.iterdir()) + ["provenance.json"],
    })
    log.info("done: %d samples (%d failed), %d genes tested, %d significant in at least one sample, "
             "consensus top %s, %.1fs, peak RSS %s MB -> %s", n_ok, len(statuses) - n_ok, len(cons),
             summary["n_genes_significant_any"], ", ".join(cons["gene"].head(5)), time.monotonic() - t0,
             peak_rss_mb(), output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
