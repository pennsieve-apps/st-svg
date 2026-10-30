# st-svg

Which genes have a spatial pattern? This app finds **spatially variable genes (SVGs)** in
spatial transcriptomics data: genes whose expression in one spot or cell is more like that of
its neighbours than chance would give. In a tissue section these are the genes that mark
structures (immunoglobulins in lymphoid aggregates and plasma-cell rich stroma, surfactant genes
in alveolar lung, keratins in tumour nests), and they are the usual starting point for describing
a tissue's architecture, choosing genes for spatial domain detection, or comparing conditions.

For each sample it computes **Moran's I** (Moran 1950) of every tested gene on the spatial
neighbour graph, with squidpy (Palla et al. 2022): I near 0 means no spatial structure, I near 1
means neighbours have almost the same expression. Optionally also **Geary's C** (Geary 1954),
which weights local differences more. P-values are analytic (fast) or from permutations, corrected
for the number of genes (Benjamini–Hochberg), and genes are ranked. A **consensus** over samples
then says in how many samples each gene is significant and how strong its autocorrelation is on
average. SpatialDE (Svensson et al. 2018) introduced SVG detection for spatial transcriptomics with
Gaussian-process models; Moran's I is the classical, much faster autocorrelation statistic
(squidpy's default for SVGs), which scales to whole Xenium slides in seconds.

It is a step of the **Spatial Domain Mapping & Condition Comparison** workflow:
`st-prep → st-qc → st-domains → {st-svg, st-compare} → st-report`, and can run after `st-qc` in
any of the three workflows. It runs once per workflow, so it **writes at the root of its single
output port**.

## Method

Per sample, streaming (one `.h5ad` in memory at a time; only `X`, `obs`/`var` names,
`obsm["spatial"]`, `uns["pennsieve"]` and the lowres image are read, not `layers["counts"]`):

1. **Expression.** `X` as st-qc wrote it: log-normalized (`normalize_total`, `log1p`). If the
   file has no `uns["pennsieve"]["normalized"]` and X holds integer counts (st-prep output), it is
   normalized here (`normalize_total` to the median, `log1p`) and provenance says so.
2. **Genes** (`GENES`): `hvg` = `var["highly_variable"]` from st-qc (all genes of a small panel
   such as Xenium's 313; when the column is absent, scanpy `highly_variable_genes(flavor="seurat",
   n_top_genes=2000)` here); `all`; or a comma-separated list. Genes detected in fewer than 3
   observations, or constant, are not tested.
3. **Spatial graph** (conventions.md, "Spatial neighbours"): Visium (`grid: hex`) uses squidpy's
   grid neighbours, 6 per spot and `N_RINGS` rings; other platforms the Delaunay triangulation or
   the k nearest neighbours (`NEIGHBOURS`), with edges longer than `MAX_NEIGHBOUR_DISTANCE` ×
   `spot_diameter` dropped so sparse regions don't connect across gaps. Binary, symmetric, then
   row-normalized by squidpy (`transformation=True`).
4. **Statistic.** `squidpy.gr.spatial_autocorr(mode="moran")`, and `mode="geary"` with
   `MODE=both`.
5. **P-values** (`N_PERMS`). `0` (default): **analytic**, the normal approximation of Moran's I
   under the normality assumption (Cliff & Ord variance as computed by squidpy, the same as
   PySAL/esda). It is accurate for the hundreds to hundreds of thousands of observations of a
   spatial sample and costs nothing. `N > 0`: permutation p-values `(k + 1)/(N + 1)` from `N`
   random relabellings (seeded by `SEED`). The test is **one-sided, for positive
   autocorrelation** (I above its expectation −1/(n−1); C below 1): squidpy reports the tail the
   statistic falls in, which would make strongly *anti*-correlated (checkerboard) genes
   "significant"; this app converts both kinds of p-value to the positive tail. Analytic p-values
   use the normal survival function, so they don't underflow to 0 until below ~1e-300.
6. **FDR.** Benjamini–Hochberg (1995) over the genes tested in the sample; `significant` =
   `pval_adj < FDR` and I > −1/(n−1). `rank` 1 = highest Moran's I in the sample.
7. **Consensus.** Per gene over the samples it was tested in: `n_samples`, `n_significant`,
   `mean_morans_i`, `median_rank`; ranked by `n_significant`, then `mean_morans_i`.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `spatial` | input (subfolder) | `<sample_id>.h5ad` | st-qc output (AnnData contract). Other files are ignored. |
| `metadata` | input (subfolder, optional) | `samples.csv` | st-prep's sample sheet; gives the sample order and `condition` (copied into summary.json). |
| `svg` | output | everything below, at the root | The only output, so it is the whole of `OUTPUT_DIR`. |

Without links the same files are found by name anywhere under `INPUT_DIR` (symlinks followed).

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `GENES` | string | `hvg` | `hvg`, `all`, or a comma-separated list of gene symbols (unknown symbols are logged). |
| `MODE` | string | `moran` | `moran`, or `both` (adds Geary's C with its own p-values). |
| `N_PERMS` | integer | `0` | `0` = analytic p-values. `N` = permutation p-values; the smallest is 1/(N+1), so with G genes and FDR q at least G/(q·(N+1)) genes must reach it before BH rejects anything (2000 HVGs at q = 0.05 need N ≳ 40 000 for a single gene). A warning is logged when N is too small. |
| `N_RINGS` | integer | `1` | Visium hex grid: rings of neighbours (1 = the 6 adjacent spots, 2 = 18). |
| `NEIGHBOURS` | string | `delaunay` | Non-grid data: `delaunay` (cells) or an integer k for k nearest neighbours (6 for Visium HD / Slide-seq bins). |
| `MAX_NEIGHBOUR_DISTANCE` | number | `3.0` | Non-grid data: edges longer than this × `spot_diameter` are dropped; `0` = no cap. |
| `FDR` | number | `0.05` | BH-adjusted p-value threshold for `significant`. |
| `TOP_N` | integer | `12` | Genes per sample in the spatial figure (at most 24 panels) and in summary.json; genes in the consensus figure. |
| `SEED` | integer | `42` | Permutation seed (`N_PERMS > 0`). |
| `THREADS` | integer | `0` | numba/BLAS threads; `0` = from the cgroup CPU quota. |

Bad parameters exit 2; no input or every sample failed exit 1; a sample that fails is logged and
skipped.

## Outputs (`OUTPUT_DIR/`)

| File | Contents |
|---|---|
| `svg.csv` | One row per sample × tested gene, sorted by sample then rank: `sample_id, gene, morans_i, pval, pval_adj, significant, rank`; with `MODE=both` also `gearys_c, gearys_pval, gearys_pval_adj`. |
| `svg_consensus.csv` | `gene, n_samples` (samples it was tested in), `n_significant, mean_morans_i, median_rank, rank`. |
| `summary.json` | `n_samples`, `n_failed`, `statistic`, `pvalues`, `fdr`, `n_genes_tested` (union), `n_genes_significant_any`, `n_genes_significant_all`, `per_sample` (`n_obs, n_genes_tested, n_significant, fraction_significant, median_morans_i, condition, top_genes, graph_mean_degree`), `consensus_top` (TOP_N genes), `parameters`, `figures`. |
| `<sample_id>_svg_top.png` | The sample's TOP_N genes by Moran's I, one panel each, log-normalized expression (1st–99th percentile colour range) on the lowres H&E in grey when there is one; I and q in the title, "(n.s.)" when not significant. |
| `morans_distribution.png` | Moran's I of every tested gene, one row per sample; significant genes dark, with the count significant / tested. |
| `consensus_top.png` | The TOP_N consensus genes × samples, Moran's I as colour, a dot where significant, grey where not tested. |
| `provenance.json` | App and tool versions, parameters, method, inputs seen, per-sample status (n_obs, normalization, graph: mode, mean degree, isolated observations, radius cap; genes selected/not tested; timings), peak RSS. |

## Results on the demo data

**TLS Visium** (st-qc output of the 8 kidney/lung tumour sections, 2 082–4 431 spots, 2 000 HVGs
each; defaults): significant SVGs per sample KC1 436, KC2 237, KC3 1 568, LC1 304, LC2 1 803,
LC3 933, LC4 493, LC5 1 613 (FDR 0.05). Top genes: immunoglobulins (IGHG1, IGHG2, IGHA1, IGKC,
IGHM: plasma cells and lymphoid aggregates), surfactant genes SFTPB/SFTPC in lung (LC2 I = 0.80,
0.78), S100A2, PRB3/PRB4 (submucosal glands, LC4), UBC and TGFBI in KC3. 4 512 of 9 260 tested
genes are significant in at least one sample; 6 in all eight (IGHG1, IGHG2, IGHA1, IGHM, IGFBP5,
TACSTD2), the consensus top four being immunoglobulins — the B-/plasma-cell programme of the
tertiary lymphoid structures this dataset was collected for.

**Xenium breast** (st-qc output, 313-gene panel): full replicates (160 674 and 116 801 cells):
312 and 311 of 313 genes significant, top KRT7 (I = 0.74 / 0.76), SCD, EPCAM, FOXA1, TACSTD2,
CEACAM6 (luminal tumour epithelium). With 20 000 cells subsampled per replicate (`xenium-20k`):
295 and 298 significant, the same top genes (KRT7 0.68 / 0.71); Moran's I is lower there because a
random 1-in-8 subsample of cells leaves neighbours further apart (mean Delaunay degree 3.0–3.6
under the 40 µm cap, ~1 300 cells isolated, against 5.7–5.8 and 250–450 on the full slides).

**Human Lymph Node** Visium (3 928 spots): IGKC, IGHG4, IGHG1, CCL21 (T-cell zone), FDCSP and
CR2 (follicular dendritic cells / germinal centres); 1 406 of 2 000 HVGs significant.

On the synthetic fixture the four planted genes rank 1–4 with I = 0.51–0.72 and q < 1e-80 in
every sample where they are patterned; at most one of 200 noise genes per sample passes FDR 0.05.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 3008 MiB, 15 min timeout. The image is 2.2 GB
(squidpy brings spatialdata, dask, geopandas, scikit-image), under Lambda's 10 GB limit, and peak
memory stays well under 3008 MiB on all the demo data, so it qualifies as a light step.

Measured in the linux/amd64 image under emulation on an M-series Mac (Docker `--memory 3g
--cpus 1 --network none`, defaults):

| Data | Observations × genes tested | Wall time | Peak RSS | cgroup `memory.peak` |
|---|---|---|---|---|
| TLS Visium, 8 samples | 2 082–4 431 spots × 2 000 | 39 s | 1.34 GB | 1.87 GB (0.68 GB page cache) |
| Xenium breast, 2 full slides | 160 674 + 116 801 cells × 313 | 43 s | 1.34 GB | 1.50 GB |
| Xenium breast, 20 k cells each | 19 189 + 19 669 × 313 | 19 s | 0.85 GB | 0.76 GB |
| TLS Visium, `MODE=both N_PERMS=999` | as above | see below | | |

About 15 s of each run is numba compiling scanpy's kernels on the first sample (emulated;
cached in `/tmp/numba` for later samples), and ~0.8 GB of the RSS is the imported stack itself.
Memory grows with one sample's `X` (sparse) and graph; a Visium HD 8 µm sample (~500 k bins ×
2 000 HVGs) should still fit in a few GiB, but `GENES=all` on such a sample, or permutations on
hundreds of thousands of cells, can exceed Lambda's 15 minutes: use the standard compute type.

Permutations cost roughly N_PERMS × the analytic test per mode; see the row above.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python    # unit tests on the synthetic fixture, no Docker
make run PYTHON=.venv/bin/python     # Docker over data/fixture/linked -> data/output/
```

By hand:

```sh
docker build --platform linux/amd64 -t st-svg .
docker run --rm --platform linux/amd64 --network none --memory 3g \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "$PWD/data/fixture/linked:/data/input:ro" -v "$PWD/data/output:/data/output" st-svg
```

The same image runs as a Lambda function: with `AWS_LAMBDA_RUNTIME_API` set,
`processor/entry.py` hands over to the runtime interface client with `processor.handler.handler`,
which takes `inputDir`, `outputDir` and the parameters from the event.

## Limitations

- **Significance is not effect size.** With thousands of spots the analytic test calls genes
  with I = 0.02–0.03 significant (half or more of the HVGs in several TLS samples). Read
  `morans_i` and `rank`, not only `significant`; the consensus ranks by how many samples agree
  and then by mean I.
- Moran's I measures **global** autocorrelation: one coherent pattern across the section. A gene
  patterned in a small region of a large section gets a small I; local statistics (LISA,
  hotspot analysis) are not computed. It is also blind to the *shape* of the pattern (gradient,
  spots, stripes all score alike) — look at the figures or use st-domains.
- **Composition drives most SVGs.** In spot data a gene is spatially variable mostly because
  the cell types expressing it are arranged in space (immunoglobulins follow plasma cells). The
  app doesn't separate cell-type composition from within-cell-type spatial regulation;
  deconvolution (cell2location) or cell-level data does that better.
- Library-size normalization can itself induce weak spatial structure in unpatterned genes where
  total counts vary across the tissue; such genes come with small I.
- Moran's I depends on the graph: more rings (`N_RINGS`), k, or a sparser subsample change its
  value. Compare I between samples only when they were built the same way (same platform and
  parameters) and at similar density.
- Permutation p-values are bounded below by 1/(N_PERMS+1); with the default 2 000 HVGs only very
  large N_PERMS can pass BH. Analytic p-values are the recommended default.
- Gene sets differ between samples when HVGs are chosen per sample; `n_samples` in the consensus
  says how many samples a gene was tested in.

## Licences

Apache-2.0 for this app (see `LICENSE`). Bundled components and their licences: `THIRD_PARTY.md`;
none is GPL.

## Citations

Moran's I (the statistic):

> Moran PAP. **Notes on continuous stochastic phenomena.** Biometrika 37(1-2):17-23 (1950).
> https://doi.org/10.1093/biomet/37.1-2.17 (not in PubMed; checked on Crossref).

Geary's C (`MODE=both`):

> Geary RC. **The Contiguity Ratio and Statistical Mapping.** The Incorporated Statistician
> 5(3):115 (1954). https://doi.org/10.2307/2986645 (not in PubMed; checked on Crossref).

Spatially variable genes in spatial transcriptomics:

> Svensson V, Teichmann SA, Stegle O. **SpatialDE: identification of spatially variable genes.**
> Nat Methods 15(5):343-346 (2018). https://doi.org/10.1038/nmeth.4636, PMID 29553579.

Squidpy (graph and statistics):

> Palla G, Spitzer H, Klein M, Fischer D, Schaar AC, Kuemmerle LB, Rybakov S, Ibarra IL,
> Holmberg O, Virshup I, Lotfollahi M, Richter S, Theis FJ. **Squidpy: a scalable framework for
> spatial omics analysis.** Nat Methods 19(2):171-178 (2022).
> https://doi.org/10.1038/s41592-021-01358-2, PMID 35102346.

False discovery rate:

> Benjamini Y, Hochberg Y. **Controlling the False Discovery Rate: A Practical and Powerful
> Approach to Multiple Testing.** J R Stat Soc Series B 57(1):289-300 (1995).
> https://doi.org/10.1111/j.2517-6161.1995.tb02031.x (not in PubMed; checked on Crossref).

Scanpy and anndata:

> Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data
> analysis.** Genome Biol 19(1):15 (2018). https://doi.org/10.1186/s13059-017-1382-0, PMID 29409532.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** JOSS 9(101):4371 (2024). https://doi.org/10.21105/joss.04371 (not in PubMed).

Demo data:

> Dawo S, Nonchev K, Silina K. **10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung
> (5) Cancer with Tertiary Lymphoid Structures.** Zenodo (2025).
> https://doi.org/10.5281/zenodo.14620362 (CC BY 4.0; checked on DataCite).

> Janesick A, Shelansky R, Gottscho AD, et al. (16 authors incl. 10x Development Teams). **High
> resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in
> situ analysis.** Nat Commun 14(1):8353 (2023). https://doi.org/10.1038/s41467-023-43458-x,
> PMID 38114474. Data: 10x Genomics, *FFPE Human Breast Cancer* (Xenium) and *Human Lymph Node*
> (Visium) dataset pages, CC BY 4.0.

Every DOI was checked on Crossref (DataCite for Zenodo) and every PMID on PubMed (E-utilities
`esearch` by DOI, then `esummary`: title, authors, year, journal, volume and pages match).
`CITATION.cff` lists every author as Crossref/DataCite give them.
