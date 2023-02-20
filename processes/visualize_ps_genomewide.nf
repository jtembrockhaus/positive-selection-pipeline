#!/usr/bin/env nextflow

/*
 * Create a list of positions under positive selection 
 * for a single gene.
 */

process VISUALIZE_PS_GENOMEWIDE {

    conda "${projectDir}/envs/visualize_ps_genomewide.yaml"

    input:
    val fel_sites_under_ps_ch
    val total_sites_ch
    val gene_lengths_ch

    output:
    path "fel_genomewide_positive_selection.html", emit: genomewide_ps
    
    script:
    """
    python ${projectDir}/scripts/visualize_ps_genomewide.py --sites_ps ${fel_sites_under_ps_ch} --sites_total ${total_sites_ch} --gene_lengths ${projectDir}/${gene_lengths_ch}
    """
}