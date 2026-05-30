#=========================
#imports
#=========================
import sys
import anndata
import pytest
import scanpy as sc
import numpy as np

sys.path.append("../source")
from normalize import normalize


#=========================
#func
#=========================

@pytest.fixture
def sample_matrix():
    return np.array([
    [1, 0, 5],
    [0, 0, 0],
    [2, 3, 1],
])

@pytest.fixture
def sample_adata(sample_matrix):
    adata = anndata.AnnData(X=sample_matrix)
    return adata

def test_normalize(sample_adata):
    result = normalize(sample_adata, target_sum=10)
    assert (result.X >= 0).all()

    