#!/usr/bin/env nextflow

/*
 * Analyze selective pressure using the Fast Unconstrained Bayesian
 * AppRoximation (FUBAR) method implemented in HyPhy.
 */

process FUBAR_ANALYSIS {

    conda "${projectDir}/envs/positive_selection_analysis.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch
    path newick_tree_ch

    output:
    path "${gene}_fubar_results.json", emit: fubar_results_ch
    
    script:
    """
    HYPHYMPI CPU=${task.cpus} fubar --grid 40 --alignment ${nuc_msa_filtered_ch} --tree ${newick_tree_ch} --output ${gene}_fubar_results.json
    """
}