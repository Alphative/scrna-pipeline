import scanpy as sc
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



MIN_GENES = 200
MAX_GENES = 2500
MAX_PCT_MT = 5.0


data_path = "../data/raw/"
adata = sc.read_10x_mtx(data_path, var_names = "gene_symbols", cache = True)

logger.info(f"Before filtering: {adata.n_obs} ")


adata.var["mt"] = adata.var_names.str.startswith("MT-")

cell_qc, gene_qc = sc.pp.calculate_qc_metrics(adata, qc_vars = ["mt"])


mask = (cell_qc["n_genes_by_counts"] > MIN_GENES) & (cell_qc["n_genes_by_counts"] < MAX_GENES) & (cell_qc["pct_counts_mt"] < MAX_PCT_MT)

adata_filtered = adata[mask]

logger.info(f"After filtering: {adata_filtered.n_obs} cells")
logger.info(f"Removed: {adata.n_obs - adata_filtered.n_obs} cells")


adata_filtered.write(filename = "../data/filtered_qc/hg_19_qc_filtered.h5ad", compression = "gzip")
