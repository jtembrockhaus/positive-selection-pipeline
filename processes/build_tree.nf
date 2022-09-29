#!/usr/bin/env nextflow

/*
 * Build a Newick tree by applying canonical neighbour-joining
 * using the software RapidNJ.
 */

process BUILD_TREE {

    conda "${projectDir}/envs/build_tree.yaml"

    input:
    val gene
    path nuc_msa_filtered_ch

    output:
    path "${gene}_newick.tree", emit: newick_tree_ch
    
    
    script:
    STO="${gene}_nuc_msa_filtered.sto"
    TREE="${gene}_newick.tree"
    """
    seqmagick convert ${nuc_msa_filtered_ch} ${STO}
    rapidnj ${STO} -i sth > ${TREE}
    sed -i "s/'//g" ${TREE}
    """
}