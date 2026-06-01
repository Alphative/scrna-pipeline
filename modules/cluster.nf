process CLUSTER {
    publishDir "${params.output_cluster}", mode: 'copy'
    
    input: path normalized_data

    output: 
    path "clustered.h5ad"
    path "metadata_clustered.json"

    script:
    """
    python /scrna_seq/source/cluster.py \
    --input_path $normalized_data \
    --output_path clustered.h5ad \
    --metadata_path metadata_clustered.json \
    --n_neighbors ${params.n_neighbors} \
    --n_comps ${params.n_comps} \
    --resolution ${params.resolution}
    """
}