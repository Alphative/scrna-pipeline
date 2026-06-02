#=========================
#imports
#=========================
import logging
import scanpy as sc
import argparse
import json
import datetime
#=========================
#logger cfg
#=========================
logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)
#=========================
#func
#=========================
def load_data(data_path):
    """Load normalized AnnData object from h5ad file."""
    adata = sc.read_h5ad(data_path)
    logger.info(f"Loaded successfully {adata.n_obs} cells")
    return adata

def run_pca(adata, n_comps=50):
    """Reduce dimensionality using PCA."""
    sc.tl.pca(adata, n_comps=n_comps)
    logger.info(f"Successfully Compressed {n_comps}")
    return adata

def run_neighbors(adata, n_neighbors=15):
    """Build k-nearest neighbor graph for clustering and UMAP."""
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    logger.info(f"Neighbors found: {n_neighbors}")
    return adata

def run_clustering(adata, resolution=0.5):
    """Cluster cells using Leiden algorithm."""
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
    logger.info(f"Clusters created: {adata.obs['leiden'].nunique()}")
    return adata

def run_umap(adata):
    """Compute UMAP embedding for visualization."""
    sc.tl.umap(adata)
    logger.info(f"UMAP created {adata.obsm['X_umap'].shape}")
    return adata

def save_results(adata, output_path):
    """Save clustered AnnData object to h5ad file."""
    adata.write(filename=output_path, compression="gzip")
    logger.info(f"File successfully saved to {output_path}")

def save_metadata(n_cells, n_clusters, metadata_path, parameters):
    """Save clustering metadata and parameters to JSON file."""
    metadata = {
        "timestamp": datetime.datetime.now().isoformat(),
        "cells": n_cells,
        "clusters": n_clusters,
        "parameters": parameters,
    }
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)
    logger.info(f"Metadata successfully saved to {metadata_path}")


#=========================
#main
#=========================
def main():
    #argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", type=str, required = True)
    parser.add_argument("--output_path", type=str, required = True)
    parser.add_argument("--metadata_path", type=str, required = True)
    parser.add_argument("--n_neighbors", type=int, default = 15)
    parser.add_argument("--n_comps", type=int, default=50)
    parser.add_argument("--resolution", type=float, default=0.5)
    args = parser.parse_args()

    pipeline_params = {
        "n_neighbors": args.n_neighbors,
        "n_comps": args.n_comps,
        "resolution": args.resolution,
    }

    #loading data
    adata = load_data(data_path=args.input_path)
    #PCA
    adata = run_pca(adata, n_comps=args.n_comps)
    #neighbors
    adata = run_neighbors(adata, n_neighbors=args.n_neighbors)
    #clustering
    adata = run_clustering(adata, resolution=args.resolution)
    #UMAP
    adata = run_umap(adata)
    #saving
    save_results(adata, args.output_path)
    save_metadata(
        n_cells=int(adata.n_obs),
        n_clusters=int(adata.obs["leiden"].nunique()),
        metadata_path=args.metadata_path,
        parameters=pipeline_params
    )

if __name__ == "__main__":
    main()