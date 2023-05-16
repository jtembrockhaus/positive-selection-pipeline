#!/usr/bin/env nextflow

/*
 * Create a codon-aware Multiple Sequence Alignment (MSA) for
 * the compressed protein sequences using MAFFT. Additionally create 
 * a position table where MSA positions are mapped to the respective 
 * reference positions.
 */

process CREATE_PROTEIN_MSA {

    conda "${projectDir}/envs/create_protein_msa.yaml"

    input:
    val gene
    path protein_seqs_compressed_ch

    output:
    path "${gene}_protein_msa.fas", emit: protein_msa_ch
    path "${gene}_position_map_table.tsv", emit: position_map_table

    script:
    REFERENCE="${projectDir}/data/static/reference_genes/reference.${gene}_protein.fas"
    """
    mafft --auto --thread ${task.cpus} --add ${protein_seqs_compressed_ch} ${REFERENCE} >| ${gene}_protein_msa.fas.tmp
    python ${projectDir}/scripts/remove_reference_from_msa.py -i ${gene}_protein_msa.fas.tmp -r ${REFERENCE} -o ${gene}_protein_msa.fas -s ${gene}_mapped_reference.fas
    python ${projectDir}/scripts/create_position_map_table.py -r ${gene}_mapped_reference.fas -p ${gene}_position_map_table.tsv
    rm ${gene}_protein_msa.fas.tmp
    """
    // DIE RICHTIGEN NEUEN OUPUTS IN main.NF bearbeiten und testen!!!
}