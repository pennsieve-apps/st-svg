# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the pinned versions in `requirements.txt`; their licences
apply to the image. **Nothing in the image is GPL** (no igraph / leidenalg: this app does no
clustering), so the image as a whole can be redistributed under these permissive licences. A
few components are MPL-2.0 (file-level copyleft: changes to those files would have to be
shared; we ship them unmodified) and one data package is CC BY 4.0 (colorcet's colour maps).

| Component | What | Licence | Source |
|---|---|---|---|
| squidpy 1.8.4 | Spatial neighbour graph (`spatial_neighbors_grid` / `_delaunay` / `_knn`), `spatial_autocorr` (Moran's I, Geary's C). | BSD-3-Clause | https://github.com/scverse/squidpy |
| scanpy 1.12.5 | `morans_i` / `gearys_c` kernels used by squidpy; normalization and HVG fallback. | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata 0.13.4 | Reading `.h5ad` elements. | BSD-3-Clause | https://github.com/scverse/anndata |
| NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.0 | Arrays, sparse matrices, normal tail probabilities, tables. | BSD-3-Clause | numpy.org, scipy.org, pandas.pydata.org |
| statsmodels 0.15.0 | Used by squidpy (multiple-testing helpers). | BSD-3-Clause | https://www.statsmodels.org |
| h5py 3.16.0 | HDF5 access. | BSD-3-Clause | https://www.h5py.org |
| numba 0.68.1, llvmlite 0.50.0 | JIT kernels of scanpy/squidpy. | BSD-2-Clause; llvmlite BSD-2-Clause AND Apache-2.0 WITH LLVM-exception | https://numba.pydata.org |
| Matplotlib 3.11.2 | PNG figures. | Matplotlib licence (PSF-based, BSD-compatible) | https://matplotlib.org |
| spatialdata 0.8.0, spatialdata-plot 0.4.3, dask 2026.8.0, distributed, xarray, zarr, ome-zarr, geopandas, shapely, pyproj, pyogrio, pyarrow, scikit-image, scikit-learn, dask-image, datashader, centrosome, cp_measure, mahotas, omnipath, ... | Dependencies of squidpy, installed but not called by this app (versions in `requirements.txt`). | BSD-2/3-Clause, MIT, Apache-2.0 (pyarrow, xarray, multiscale_spatial_image, docrep, ...) | PyPI |
| certifi, fast-array-utils, legacy-api-wrap, session-info2, tqdm | Dependencies (versions in `requirements.txt`). | MPL-2.0 (tqdm: MPL-2.0 AND MIT) | PyPI |
| colorcet 3.2.1 | Dependency of datashader. | CC-BY-4.0 | https://github.com/holoviz/colorcet |
| Other transitive dependencies (aiohttp, attrs, click, cloudpickle, fsspec, Jinja2, joblib, networkx, packaging, pillow, pydantic, PyYAML, requests, rich, ...) | Versions in `requirements.txt`. | MIT / BSD / Apache-2.0 / PSF-2.0 / MIT-CMU (Pillow) | PyPI |
| awslambdaric 4.2.3 (+ simplejson 4.2.0, snapshot-restore-py 1.0.0) | AWS Lambda runtime interface client, for the Lambda runtime. | Apache-2.0 (simplejson: MIT OR AFL-2.1) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12.11-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

Licences were read from each installed package's metadata in the built image
(`importlib.metadata`), and from the LICENSE file where the metadata had none (cp_measure: MIT).

## Not bundled

- **SpatialDE** (Svensson et al. 2018) is cited as the reference method for spatially variable
  genes; it is not used or bundled (this app uses Moran's I, a faster autocorrelation statistic).
- **Demo data** (not in the image or the repository): the TLS Visium dataset (Dawo, Nonchev,
  Silina; Zenodo 2025, doi:10.5281/zenodo.14620362, CC BY 4.0) and the 10x Genomics Xenium FFPE
  Human Breast Cancer replicates (CC BY 4.0).
