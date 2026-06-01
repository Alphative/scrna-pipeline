    process NORMALIZE {
        publishDir "${params.output_norm}", mode: 'copy'
        input: path filtered_data

        output: 
        path "normalized.h5ad"
        path "metadata_norm.json"

        script:
        """
        python normalize.py \
        --input_path $filtered_data \
        --output_path normalized.h5ad \
        --metadata_path metadata_norm.json \
        --target_sum ${params.target_sum}
        """
    }