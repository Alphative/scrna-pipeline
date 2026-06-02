# =========================
# imports
# =========================
import argparse
import datetime
import json
import logging
import scanpy as sc

# =========================
# logging setup
# =========================
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =========================
# func
# =========================
def load_data(data_path):
    """Load 10x MTX matrix from directory."""
    adata = sc.read_10x_mtx(data_path, var_names="gene_symbols", cache=True)
    logger.info(f"Before filtering: {adata.n_obs} ")
    return adata


def compute_qc(adata, mt_prefix):
    """Compute QC metrics: mitochondrial, ribosomal, and MALAT1 gene tracking."""
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["malat1"] = adata.var_names == "MALAT1"
    logger.info(f"Tracking mitochondrial genes: {adata.var['mt'].sum()}")
    logger.info(f"Tracking ribosomal genes: {adata.var['ribo'].sum()}")
    logger.info(f"Tracking MALAT1: {adata.var['malat1'].sum()}")
    cell_qc, gene_qc = sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "malat1"], percent_top=None
    )

    return cell_qc, gene_qc


def filter_genes(adata):
    """Remove genes not expressed in any cell."""
    logger.info(f"Genes Before: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=1)
    logger.info(f"Genes After:{adata.n_vars}")
    return adata


def filter_cells(adata, cell_qc, min_genes, max_genes, max_pct_mt):
    """Filter cells by gene count and mitochondrial RNA percentage."""
    mask = (
        (cell_qc["n_genes_by_counts"] > min_genes)
        & (cell_qc["n_genes_by_counts"] < max_genes)
        & (cell_qc["pct_counts_mt"] < max_pct_mt)
    )
    adata_filtered = adata[mask].copy()
    logger.info(f"After filtering: {adata_filtered.n_obs} cells")
    logger.info(f"Removed: {adata.n_obs - adata_filtered.n_obs} cells")
    return adata_filtered


def save_results(adata_filtered, output_path):
    """Save filtered AnnData object to h5ad file."""
    adata_filtered.write(filename=output_path, compression="gzip")
    logger.info(f"File successfully saved to {output_path}")


def save_metadata(input_cells, output_cells, output_path, parameters):
    """Save pipeline run metadata and parameters to JSON file."""
    metadata = {
        "timestamp": datetime.datetime.now().isoformat(),
        "input_cells": input_cells,
        "output_cells": output_cells,
        "removed_cells": input_cells - output_cells,
        "parameters": parameters,
    }
    with open(output_path, "w") as f:
        json.dump(metadata, f, indent=4)


# =========================
# MAIN
# =========================
def main():
    # agrparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_genes", type=int, default=200)
    parser.add_argument("--max_genes", type=int, default=2500)
    parser.add_argument("--max_pct_mt", type=float, default=5.0)
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--output_path", type=str, required=True)
    parser.add_argument("--metadata_path", type=str, required=True)
    parser.add_argument("--mt_prefix", type=str, default="MT-")
    args = parser.parse_args()

    pipeline_params = {
        "min_genes": args.min_genes,
        "max_genes": args.max_genes,
        "max_pct_mt": args.max_pct_mt,
        "mt_prefix": args.mt_prefix,
    }

    # loading data
    adata = load_data(data_path=args.data_path)
    input_cells_count = int(adata.n_obs)

    # filtering dead genes
    adata = filter_genes(adata)

    # compute qc
    cell_qc, gene_qc = compute_qc(adata, mt_prefix=args.mt_prefix)

    # filtering
    adata_filtered = filter_cells(
        adata,
        cell_qc,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_pct_mt=args.max_pct_mt,
    )

    # saving
    save_results(adata_filtered, args.output_path)
    save_metadata(
        input_cells=input_cells_count,
        output_cells=int(adata_filtered.n_obs),
        output_path=args.metadata_path,
        parameters=pipeline_params,
    )


if __name__ == "__main__":
    main()