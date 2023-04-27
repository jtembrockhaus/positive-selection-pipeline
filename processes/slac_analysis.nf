#!/usr/bin/env nextflow

/*
 * Analyze selective pressure using the Single-Likelihood Ancestor
 * Counting (SLAC) method implemented in HyPhy.
 */

process SLAC_ANALYSIS {

    conda "${projectDir}/envs/positive_selection_analysis.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch
    path newick_tree_ch

    output:
    path "${gene}_slac_results.json", emit: slac_results_ch
    
    script:
    """
    HYPHYMPI CPU=${task.cpus} slac --alignment ${nuc_msa_filtered_ch} --tree ${newick_tree_ch} --branches All --samples 0 --output ${gene}_slac_results.json
    """
}