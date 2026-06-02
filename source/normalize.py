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
    """Load filtered AnnData object from h5ad file"""
    adata = sc.read_h5ad(data_path)
    logger.info(f"Loaded successfully {adata.n_obs} cells")
    return adata

def normalize(adata, target_sum=1e4):
    """Normalize counts per cell to target_sum and apply log1p transformation"""
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    logger.info(f"Normalized {adata.n_obs} cells successfully")
    return adata

def save_results(adata, output_path):
    """Save normalized AnnData object to h5ad file."""
    adata.write(filename=output_path, compression="gzip")
    logger.info(f"File successfully saved to {output_path}")


def save_metadata(adata, metadata_path, target_sum):
    """Save normalization metadata and parameters to JSON file."""
    metadata = {"timestamp": datetime.datetime.now().isoformat(),
                "cells": adata.n_obs, 
                "parameters": {"target_sum": target_sum}}
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4)


#=========================
#main
#=========================
def main():
    #argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", type=str, required = True)
    parser.add_argument("--output_path", type=str, required = True)
    parser.add_argument("--metadata_path", type=str, required = True)
    parser.add_argument("--target_sum", type=float, default = 1e4)
    args = parser.parse_args()
    #loading data
    adata = load_data(data_path=args.input_path)
    #normalizing
    adata = normalize(adata, target_sum=args.target_sum)
    #saving
    save_results(adata, args.output_path)
    save_metadata(adata, args.metadata_path, target_sum=args.target_sum)
if __name__ == "__main__":
    main()
