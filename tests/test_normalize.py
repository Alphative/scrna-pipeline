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
], dtype=np.float32)

@pytest.fixture
def sample_adata(sample_matrix):
    adata = anndata.AnnData(X=sample_matrix)
    return adata

def test_normalize(sample_adata):
    target_sum = 10
    result = normalize(sample_adata, target_sum=10)
    assert result.X.shape == (3, 3)
    assert (result.X >= 0).all()
    assert result.X[1].sum() == 0.0
    
    expected_val = np.log1p((1.0 / 6.0) * target_sum)
    assert result.X[0, 0] == pytest.approx(expected_val, rel=1e-5)

    