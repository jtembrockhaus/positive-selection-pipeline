#!/usr/bin/env nextflow

/*
 * Create a figure showing the evolution of a single site over the time.
 */

process TEMPORAL_EVOLUTION_PLOT {

    conda "${projectDir}/envs/temporal_evolution_plot.yaml"

    input:
    val genes_ch
    val single_sites_ch
    val protein_msa_ch
    val metadata_ch
    val protein_duplicates_ch
    
    script:
    """
    python ${projectDir}/scripts/temporal_evolution_plot.py \
    --genes ${genes_ch} \
    --sites ${single_sites_ch} \
    --msa ${protein_msa_ch} \
    --metadata ${metadata_ch} \
    --duplicates ${protein_duplicates_ch}
    """
}