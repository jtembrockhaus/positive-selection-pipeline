#!/usr/bin/env nextflow

/*
 * Create a report file for the pipeline run.
 */

process PIPELINE_REPORT {

    conda "${projectDir}/envs/pipeline_report.yaml"

    input:
    val gene
    path copies_ch
    path newick_tree_ch
    path nuc_msa_filtered_ch

    output:
    path "pipeline_report.txt", emit: pipeline_report
    
    script:
    """
    python ${projectDir}/scripts/pipeline_report.py \
    --gene ${gene} \
    --pipeline_dir ${projectDir} \
    --duplicated_seqs ${copies_ch} \
    --tree ${newick_tree_ch} \
    --msa ${nuc_msa_filtered_ch} \
    --output_file pipeline_report.txt
    """
}