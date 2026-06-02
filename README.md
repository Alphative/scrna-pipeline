# About

This project is an automated cloud-ready pipeline for single cell RNA sequencing (scRNA-seq) data analysis, consisting of three stages: Quality Control, Normalization, and Clustering

## Requirements

- Docker
- Nextflow
- AWS CLI (optional, for cloud deployment)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Alphative/scrna-pipeline
```

2. Build the Docker image:
```bash
docker build -t scrna-pipeline -f docker/Dockerfile .
```

## Running locally

```bash
nextflow run main.nf
```

With custom parameters:
```bash
nextflow run main.nf --min_genes 200 --max_genes 2500
```

Use `-resume` to continue an interrupted run.

## Running on AWS

1. Create an S3 bucket and upload your data:
```bash
aws s3 cp data/raw/ s3://your-bucket/raw/ --recursive
```

2. Update `nextflow.config` params with your S3 paths:
```nextflow
data_path = "s3://your-bucket/raw/"
output_qc = "s3://your-bucket/results/qc/"
```

3. Run with AWS Batch profile:
```bash
nextflow run main.nf -profile awsbatch
```

## Project Structure

```
scrna_pipeline/
├── main.nf              # Main Nextflow pipeline
├── nextflow.config      # Pipeline configuration
├── modules/             # Nextflow process modules
│   ├── qc.nf
│   ├── normalize.nf
│   └── cluster.nf
├── source/              # Python scripts
│   ├── qc.py
│   ├── normalize.py
│   └── cluster.py
├── tests/               # pytest tests
├── docker/              # Dockerfile
└── data/                # Input/output data (not tracked)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_genes` | 200 | Minimum genes per cell |
| `max_genes` | 2500 | Maximum genes per cell |
| `max_pct_mt` | 5.0 | Maximum mitochondrial % |
| `mt_prefix` | MT- | Mitochondrial gene prefix |
| `target_sum` | 10000 | Normalization target |
| `n_neighbors` | 15 | Neighbors for graph |
| `n_comps` | 50 | PCA components |
| `resolution` | 0.5 | Clustering resolution |