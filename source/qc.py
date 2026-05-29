#imports
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import argparse



#logging setup

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


#agrparse
parser = argparse.ArgumentParser()
parser.add_argument("--min_genes", type = int, default=200)
parser.add_argument("--max_genes", type = int, default=2500)
parser.add_argument("--max_pct_mt", type = float, default=5.0)
parser.add_argument("--data_path", type = str, required=True)
parser.add_argument("--output_path", type = str, required=True)
args=parser.parse_args()



#func

def load_data(data_path):
    """Load 10x MTX matrix from directory."""
    adata = sc.read_10x_mtx(data_path, var_names = "gene_symbols", cache = True)
    logger.info(f"Before filtering: {adata.n_obs} ")
    return adata

def compute_qc(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    cell_qc, gene_qc = sc.pp.calculate_qc_metrics(adata, qc_vars = ["mt"])
    return cell_qc, gene_qc

def filter_cells(adata, cell_qc, min_genes, max_genes, max_pct_mt):
    mask = (cell_qc["n_genes_by_counts"] > min_genes) & (cell_qc["n_genes_by_counts"] < max_genes) & (cell_qc["pct_counts_mt"] < max_pct_mt)
    adata_filtered = adata[mask]
    logger.info(f"After filtering: {adata_filtered.n_obs} cells")
    logger.info(f"Removed: {adata.n_obs - adata_filtered.n_obs} cells")
    return adata_filtered

def save_results(adata_filtered, output_path):
    adata_filtered.write(filename = output_path, compression = "gzip")
    logger.info(f"File succesfully saved to {output_path}")

#MAIN
def main():
    adata = load_data(data_path=args.data_path)
    cell_qc, gene_qc = compute_qc(adata)
    adata_filtered = filter_cells(adata, cell_qc, min_genes=args.min_genes, max_genes=args.max_genes, max_pct_mt=args.max_pct_mt)
    save_results(adata_filtered, args.output_path)

if __name__ == "__main__" :
    main()