#!/usr/bin/env nextflow

/*
 * Create figures and tables of the selection analysis results.
 */

process VISUALIZE_RESULTS {

    conda "${projectDir}/envs/visualize_results.yaml"

    input:
    val gene
    path fel_results_ch
    //path meme_results_ch
    
    script:
    // """
    // python ${projectDir}/scripts/visualize_results.py --gene ${gene} --fel ${fel_results_ch} --meme ${meme_results_ch}
    // """
    """
    python ${projectDir}/scripts/visualize_results.py --gene ${gene} --fel ${fel_results_ch}
    """
}