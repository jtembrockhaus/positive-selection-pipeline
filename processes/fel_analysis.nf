#!/usr/bin/env nextflow

/*
 * Analyze selective pressure using the Fixed Efects Likelihood (FEL)
 * method implemented in HyPhy.
 */

process FEL_ANALYSIS {

    conda "${projectDir}/envs/positive_selection_analysis.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch
    path newick_tree_ch

    output:
    path "${gene}_fel_results.json", emit: fel_results_ch
    
    script:
    """
    hyphy CPU=${task.cpus} fel --alignment ${nuc_msa_filtered_ch} --tree ${newick_tree_ch} --branches Internal --output ${gene}_fel_results.json
    """
}