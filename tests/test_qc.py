#
#import
#

import sys
import anndata
import numpy as np
import pandas as pd
import pytest

sys.path.append("../source")
from qc import compute_qc, filter_cells, load_data

#
#cfg
#

#
#func
#
@pytest.fixture
def sample_matrix():
    return np.array([
    [1, 0, 5],
    [0, 0, 0],
    [2, 3, 1],
])

@pytest.fixture
def sample_adata(sample_matrix):
    adata = anndata.AnnData(
    X=sample_matrix,
    var=pd.DataFrame(index=["MT-ND1", "GAPDH", "MT-CO1"])
)
    return adata

def test_compute_qc(sample_adata):
    cell_qc, gene_qc = compute_qc(sample_adata)
    assert len(cell_qc) == 3
    assert "pct_counts_mt" in cell_qc.columns
    assert len(gene_qc) == 3


def test_filter_cells(sample_adata):
    cell_qc, gene_qc = compute_qc(sample_adata)
    min_genes = 1
    max_genes = 10
    max_pct_mt = 51.0
    test_filtered_cells = filter_cells(sample_adata, cell_qc, min_genes, max_genes, max_pct_mt)
    print(cell_qc[["n_genes_by_counts", "pct_counts_mt"]])   
    assert len(test_filtered_cells) == 1


#
#main
#
