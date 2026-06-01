process QC {
    publishDir "${params.output_path}", mode: 'copy'
    input: path data_path

    output:
    path "filtered.h5ad"
    path "metadata_qc.json"
    
    script:
    """
    python qc.py \
    --data_path $data_path \
    --output_path filtered.h5ad \
    --metadata_path metadata_qc.json \
    --min_genes ${params.min_genes} \
    --max_genes ${params.max_genes} \
    --max_pct_mt ${params.max_pct_mt} \
    --mt_prefix ${params.mt_prefix}
    """
}

