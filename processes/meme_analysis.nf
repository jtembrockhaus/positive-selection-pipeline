#!/usr/bin/env nextflow

/*
 * Analyze selective pressure using the Mixed Effects Model of 
 * Evolution (MEME) method implemented in HyPhy.
 */

process MEME_ANALYSIS {

    conda "${projectDir}/envs/positive_selection_analysis.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch
    path newick_tree_ch

    output:
    path "${gene}_meme_results.json", emit: meme_results_ch
    
    script:
    """
    HYPHYMPI CPU=${task.cpus} meme --alignment ${nuc_msa_filtered_ch} --tree ${newick_tree_ch} --branches Internal --output ${gene}_meme_results.json
    """
}