#!/usr/bin/env nextflow

/*
 * Create a codon-aware Multiple Sequence Alignment (MSA) for
 * the compressed protein sequences using Virulign.
 */

process CREATE_PROTEIN_MSA {

    conda "${projectDir}/envs/create_protein_msa.yaml"

    input:
    val gene
    path protein_seqs_compressed_ch

    output:
    path "${gene}_protein_msa.fas", emit: protein_msa_ch

    script:
    REFERENCE="${projectDir}/data/static/reference_genes/reference.${gene}_protein.fas"
    """
    mafft --auto --thread ${task.cpus} --add ${protein_seqs_compressed_ch} ${REFERENCE} >| ${gene}_protein_msa.fas.tmp
    python ${projectDir}/scripts/remove_reference_from_msa.py -i ${gene}_protein_msa.fas.tmp -r ${REFERENCE} -o ${gene}_protein_msa.fas
    rm ${gene}_protein_msa.fas.tmp
    """
}