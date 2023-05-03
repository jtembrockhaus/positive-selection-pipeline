#!/usr/bin/env nextflow

/*
 * Analyze selective pressure using the PRIME method implemented in HyPhy.
 */

process PRIME_ANALYSIS {

    conda "${projectDir}/envs/positive_selection_analysis.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch
    path newick_tree_ch

    output:
    path "${gene}_prime_results.json", emit: prime_results_ch
    
    script:
    """
    hyphy CPU=${task.cpus} prime --alignment ${nuc_msa_filtered_ch} --tree ${newick_tree_ch} --branches Internal --output ${gene}_prime_results.json
    """
}