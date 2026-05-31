#=========================
#imports
#=========================
import os
import sys
import anndata
import pytest
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
source_path = os.path.abspath(os.path.join(current_dir, "../source"))
if source_path not in sys.path:
    sys.path.append(source_path)
from cluster import run_pca, run_neighbors, run_clustering, run_umap

@pytest.fixture
def sample_matrix():
    return np.random.randint(0, 10, size=(30, 50))

@pytest.fixture
def sample_adata(sample_matrix):
    adata = anndata.AnnData(X=sample_matrix)
    return adata


def test_run_pca(sample_adata):
    n_comps=1
    sample_adata = run_pca(sample_adata, n_comps=n_comps)
    assert sample_adata.obsm["X_pca"].shape[1] == n_comps

def test_run_neighbors(sample_adata):
    n_neighbors=2
    sample_adata = run_pca(sample_adata, n_comps=1)
    sample_adata = run_neighbors(sample_adata, n_neighbors=n_neighbors)
    assert "connectivities" in sample_adata.obsp

def test_run_clustering(sample_adata):
    n_comps=1
    n_neighbors=2
    resolution=0.5
    sample_adata = run_pca(sample_adata, n_comps=n_comps)
    sample_adata = run_neighbors(sample_adata, n_neighbors=n_neighbors)
    sample_adata = run_clustering(sample_adata, resolution=resolution)
    assert "leiden" in sample_adata.obs

def test_run_umap(sample_adata):
    n_comps=1
    n_neighbors=2
    sample_adata = run_pca(sample_adata, n_comps=n_comps)
    sample_adata = run_neighbors(sample_adata, n_neighbors=n_neighbors)
    sample_adata = run_umap(sample_adata)
    assert "X_umap" in sample_adata.obsm