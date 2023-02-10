#!/usr/bin/env nextflow

/*
 * 
 */

process EXTRACT_EVOLUTIONARY_ANNOTATION {

    conda "${projectDir}/envs/extract_evolutionary_annotation.yaml"

    input:
    val gene
    path prime_results_ch
    path nuc_msa_merged_ch

    output:
    path "${gene}_fubar_results.json", emit: fubar_results_ch
    
    script:
    """
    python ${projectDir}/scripts/extract_evolutionary_annotation.py \\
    --region ${gene} \\
    --prime ${prime_results_ch} \\
    --coordinates ${nuc_msa_merged_ch} \\
    --offset 180 \\
    --output ${gene}_evo_annotation.json

    """
}