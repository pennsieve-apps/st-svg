"""
Figures (PNG, 150 dpi):

  <sample_id>_svg_top.png   the sample's top genes by Moran's I, one panel per gene, on the
                            lowres H&E (shown in grey so the expression colour stands out)
                            when there is one; y inverted (image convention)
  morans_distribution.png   Moran's I of every tested gene, one row per sample, significant
                            genes dark, with the count significant / tested
  consensus_top.png         the consensus top genes x samples, Moran's I as colour, a dot
                            where the gene is significant in that sample
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

DPI = 150
# One-hue sequential ramp (light -> dark blue) for magnitudes.
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#b8b7b1"
SURFACE = "#fcfcfb"
MAX_PANELS = 24


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
    return plt


def _cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("ramp", RAMP)


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def _grey(image: np.ndarray) -> np.ndarray:
    img = np.asarray(image)
    if img.dtype != np.uint8:
        img = np.clip(img * (255 if img.max() <= 1.0 else 1), 0, 255).astype(np.uint8)
    if img.ndim == 3:
        img = img[..., :3].astype(np.float32) @ np.array([0.299, 0.587, 0.114], dtype=np.float32)
    # lift towards white so the coloured spots carry the contrast
    return 255 - (255 - img.astype(np.float32)) * 0.55


def spatial_top(sample_id: str, coords: np.ndarray, values: list[np.ndarray], genes: list[str],
                stats: list[str], image: np.ndarray | None, scale: float | None, diameter: float | None,
                path: Path) -> Path:
    plt = _plt()
    n = min(len(genes), MAX_PANELS)
    ncol = min(4, max(1, n))
    nrow = max(1, math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow + 0.5), squeeze=False)
    xy = coords * (scale if (image is not None and scale) else 1.0)
    d = (diameter or 0) * (scale if (image is not None and scale) else 1.0)
    pad = max(d, 0.02 * float(np.ptp(xy[:, 0]) or 1))
    x0, x1 = xy[:, 0].min() - pad, xy[:, 0].max() + pad
    y0, y1 = xy[:, 1].min() - pad, xy[:, 1].max() + pad
    from matplotlib.ticker import MaxNLocator

    cmap = _cmap()
    used = []
    for k, ax in enumerate(axes.ravel()):
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if k >= n:
            ax.set_visible(False)
            continue
        if image is not None:
            ax.imshow(_grey(image), cmap="gray", vmin=0, vmax=255, interpolation="nearest")
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)            # y down, image convention
        ax.set_aspect("equal")
        ax.set_title(f"{genes[k]}\n{stats[k]}", fontsize=8, color=INK)
        used.append(ax)
    fig.suptitle(f"{sample_id}: top spatially variable genes (log-normalized expression, 1st-99th percentile)",
                 fontsize=10, color=INK)
    height = 3.0 * nrow + 0.5
    fig.tight_layout(w_pad=5.0, rect=(0, 0, 1 - 0.4 / ncol, 1 - 0.3 / height))
    # Marker size from the final axes geometry, so a spot is drawn at its diameter.
    used[0].apply_aspect()
    bbox = used[0].get_window_extent()
    pts_per_unit = (bbox.width * 72.0 / fig.dpi) / max(x1 - x0, 1e-9)
    size = min((max(d * pts_per_unit, 1.0)) ** 2 if d else 1.5, 80.0)
    for k, ax in enumerate(used):
        v = values[k]
        # 1st to 99th percentile, so a gene expressed everywhere still shows its pattern
        vmin, vmax = (float(np.quantile(v, 0.01)), float(np.quantile(v, 0.99))) if v.size else (0.0, 1.0)
        if vmax <= vmin:
            vmin, vmax = float(v.min()), (float(v.max()) if v.max() > v.min() else float(v.min()) + 1.0)
        order = np.argsort(v, kind="mergesort")   # high values drawn on top
        sc = ax.scatter(xy[order, 0], xy[order, 1], c=v[order], s=size, cmap=cmap, vmin=vmin, vmax=vmax,
                        linewidths=0, alpha=0.95, rasterized=True)
        cax = ax.inset_axes([1.02, 0.1, 0.035, 0.8])
        cb = fig.colorbar(sc, cax=cax)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=6, length=2)
        cb.locator = MaxNLocator(4)
        cb.update_ticks()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def morans_distribution(per_sample: dict[str, tuple[np.ndarray, np.ndarray]], path: Path) -> Path:
    """per_sample: sample_id -> (morans_i, significant)."""
    plt = _plt()
    ids = list(per_sample)
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(ids) + 1.4))
    rng = np.random.default_rng(0)
    for row, sid in enumerate(ids):
        i, sig = per_sample[sid]
        jitter = rng.uniform(-0.28, 0.28, size=i.size)
        ax.scatter(i[~sig], row + jitter[~sig], s=3, color=MUTED, linewidths=0, rasterized=True,
                   label="not significant" if row == 0 else None)
        ax.scatter(i[sig], row + jitter[sig], s=9, color=RAMP[5], linewidths=0, rasterized=True,
                   label="significant" if row == 0 else None)
        ax.text(1.01, row, f"{int(sig.sum())} / {i.size}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=7, color=INK_2)
    ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels(ids)
    ax.set_ylim(len(ids) - 0.5, -0.5)
    ax.set_xlabel("Moran's I")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, markerscale=2, fontsize=7)
    ax.set_title("Moran's I per gene and sample (right: significant / tested)\n", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def consensus_top(genes: list[str], samples: list[str], values: np.ndarray, significant: np.ndarray,
                  n_sig: list[int], path: Path) -> Path:
    """values, significant: (genes, samples); NaN where the gene was not tested in the sample."""
    plt = _plt()
    from matplotlib.colors import Normalize

    fig, ax = plt.subplots(figsize=(1.6 + 0.55 * len(samples) + 1.2, 0.28 * len(genes) + 1.8))
    cmap = _cmap().with_extremes(bad="#ebeae6")
    finite = values[np.isfinite(values)]
    vmax = float(finite.max()) if finite.size and finite.max() > 0 else 1.0
    im = ax.imshow(np.ma.masked_invalid(np.clip(values, 0, None)), cmap=cmap, norm=Normalize(0, vmax),
                   aspect="auto", interpolation="nearest")
    yy, xx = np.nonzero(significant)
    ax.scatter(xx, yy, s=6, color=SURFACE, linewidths=0)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"{g}  ({n}/{len(samples)})" for g, n in zip(genes, n_sig)])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label("Moran's I (below 0 shown as 0)")
    cb.outline.set_visible(False)
    ax.set_title("Consensus top spatially variable genes\ndot: significant in the sample; grey: not tested\n"
                 "(n/N): significant in n of N samples", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path
