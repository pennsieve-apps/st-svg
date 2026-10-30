# sample-input

This app reads the output of `st-qc` (port `spatial`: `<sample_id>.h5ad`, X log-normalized) and,
optionally, `st-prep`'s `metadata/samples.csv`. It also runs on `st-prep` output directly (raw
counts in X are normalized here, HVGs computed here).

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/
```

Three samples following the AnnData contract after st-qc: `vis_a`, `vis_b` (Visium-like hex
grids, 480 spots, a lowres "H&E" image) and `cells_a` (600 Xenium-like cells in microns), 206
genes. Planted patterns: `GRAD` (gradient), `BLOB` (disc), `STRIPE` (band), `PATCHY` (two discs;
flat in `vis_b`). Everything else is unpatterned Poisson noise; `RARE` is detected in 2
observations (not tested) and `N199` is not an HVG. Expected: the four planted genes are ranked
1-4 and significant (PATCHY not in `vis_b`), noise genes rarely so; consensus top three GRAD,
BLOB, STRIPE.

## Real data

Any spatial transcriptomics data that `st-prep` reads, through `st-prep` then `st-qc`. The demo
data of the workflows:

- **TLS Visium** (8 FFPE kidney and lung tumour sections): Dawo S, Nonchev K, Silina K. Zenodo
  2025, https://doi.org/10.5281/zenodo.14620362 (CC BY 4.0).
- **Xenium FFPE Human Breast Cancer** Rep1/Rep2 (10x Genomics, CC BY 4.0; Janesick et al. 2023).
- **Human Lymph Node** Visium (10x Genomics, CC BY 4.0).

`scripts/fetch-sample-data.sh` in the workflow repository downloads them. Then:

```sh
# st-prep -> st-qc (see those apps), then:
docker run --rm --platform linux/amd64 --network none --memory 3g \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "$WORK/tls/qc/spatial:/data/input/spatial:ro" \
  -v "$WORK/tls/prep/metadata:/data/input/metadata:ro" \
  -v "$WORK/tls/svg:/data/output" st-svg
```
