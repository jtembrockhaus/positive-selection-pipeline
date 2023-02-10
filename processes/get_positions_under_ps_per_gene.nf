#!/usr/bin/env nextflow

/*
 * Create a list of positions under positive selection 
 * for a single gene.
 */

process GET_POSITIONS_UNDER_PS_PER_GENE {

    conda "${projectDir}/envs/get_positions_under_ps_per_gene.yaml"

    input:
    val gene
    path fel_results_ch

    output:
    path "${gene}_fel_sites_under_ps.txt", emit: fel_sites_under_ps
    path "${gene}_total_sites.txt", emit: total_sites
    
    script:
    """
    python ${projectDir}/scripts/get_positions_under_ps_per_gene.py --gene ${gene} --fel ${fel_results_ch}
    """
}