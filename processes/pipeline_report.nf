#!/usr/bin/env nextflow

/*
 * Create a report file for the pipeline run.
 */

process PIPELINE_REPORT {

    conda "${projectDir}/envs/pipeline_report.yaml"

    input:
    val gene
    path fel_results_ch
    path meme_results_ch
    
    script:
    """
    python ${projectDir}/scripts/pipeline_report.py
    """
}