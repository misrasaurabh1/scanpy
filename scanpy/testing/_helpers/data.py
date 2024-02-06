"""
Functions returning copies of datasets as cheaply as possible,
i.e. without having to hit the disk or (in case of ``_pbmc3k_normalized``) recomputing normalization.
"""

from __future__ import annotations

import warnings

try:
    from functools import cache
except ImportError:  # Python < 3.9
    from functools import lru_cache

    def cache(func):
        return lru_cache(maxsize=None)(func)


from typing import TYPE_CHECKING

import scanpy as sc

if TYPE_CHECKING:
    from anndata import AnnData

# Functions returning the same objects (easy to misuse)


_pbmc3k = cache(sc.datasets.pbmc3k)
_pbmc3k_processed = cache(sc.datasets.pbmc3k_processed)
_pbmc68k_reduced = cache(sc.datasets.pbmc68k_reduced)
_krumsiek11 = cache(sc.datasets.krumsiek11)
_paul15 = cache(sc.datasets.paul15)


# Functions returning copies


def pbmc3k() -> AnnData:
    return _pbmc3k().copy()


def pbmc3k_processed() -> AnnData:
    return _pbmc3k_processed().copy()


def pbmc68k_reduced() -> AnnData:
    return _pbmc68k_reduced().copy()


def krumsiek11() -> AnnData:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "Observation names are not unique", module="anndata"
        )
        return _krumsiek11().copy()


def paul15() -> AnnData:
    return _paul15().copy()


@cache
def _pbmc3k_normalized() -> AnnData:
    pbmc = pbmc3k()
    sc.pp.filter_genes(pbmc, min_counts=1)
    sc.pp.log1p(pbmc)
    sc.pp.normalize_total(pbmc)
    sc.pp.highly_variable_genes(pbmc, inplace=True)
    return pbmc


def pbmc3k_normalized() -> AnnData:
    return _pbmc3k_normalized().copy()


def pbmc3k() -> AnnData:
    return _pbmc3k()


def filter_genes(
    data: AnnData | spmatrix | np.ndarray,
    *,
    min_counts: int | None = None,
    min_cells: int | None = None,
    max_counts: int | None = None,
    max_cells: int | None = None,
    inplace: bool = True,
) -> AnnData | tuple[np.ndarray, np.ndarray] | None:

    n_given_options = sum(
        option is not None for option in [min_cells, min_counts, max_cells, max_counts]
    )
    if n_given_options != 1:
        raise ValueError(
            "Only provide one of the optional parameters `min_counts`, "
            "`min_cells`, `max_counts`, `max_cells` per call."
        )

    adata = data if inplace else data.copy()
    X = adata.X if isinstance(data, AnnData) else data

    min_number = min_counts if min_cells is None else min_cells
    max_number = max_counts if max_cells is None else max_cells
    number_per_gene = np.sum(X > 0 if (min_cells or max_cells) else X, axis=0)

    if min_number is not None:
        gene_subset = number_per_gene >= min_number
    if max_number is not None:
        gene_subset = number_per_gene <= max_number

    # ... (rest of the function remains unchanged)


def log1p(
    data: AnnData | np.ndarray | spmatrix,
    *,
    base: Number | None = None,
    chunked: bool | None = None,
    chunk_size: int | None = None,
    layer: str | None = None,
    obsm: str | None = None,
) -> AnnData | np.ndarray | spmatrix | None:

    _check_array_function_arguments(
        chunked=chunked, chunk_size=chunk_size, layer=layer, obsm=obsm
    )
    X = data.X
    if base:
        X = np.log1p(X, out=X) / np.log(base)
    else:
        np.log1p(X, out=X)
    return data


def normalize_total(
    adata: AnnData,
    *,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: str | None = None,
    layer: str | None = None,
    layer_norm: str | None = None,
    inplace: bool = True,
) -> AnnData | dict[str, np.ndarray] | None:

    # ...(most of the function remains unchanged)

    # Update counts and data in one go
    counts_per_cell = counts_per_cell if inplace else counts_per_cell.copy()
    _normalize_data(X, counts_per_cell, target_sum, copy=not inplace)

    return adata if inplace else {"X": X, "norm_factor": counts_per_cell}
