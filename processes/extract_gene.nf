#!/usr/bin/env nextflow

/*
 * Extract the appropriate gene segment from input genomes 
 * using codon-aware alignment, filter based on N content 
 * and translate to amino-acids.
 */

process EXTRACT_GENE {

    conda "${projectDir}/envs/extract_gene.yaml"

    input:
    val sequences_ch
    val metadata_ch
    val gene
    val trim_from
    val trim_to
    val n_frac

    output:
    path "${gene}_protein.fas", emit: protein_seqs_ch
    path "${gene}_nuc.fas", emit: nuc_seqs_ch

    script:
    PREMSA="${projectDir}/ressources/hyphy-analyses/codon-msa/pre-msa.bf"
    REFERENCE="${projectDir}/data/static/reference_genes/${gene}.fas"
    TMP_FILE="${gene}"
    """
    cp ${sequences_ch} ${TMP_FILE}
    HYPHYMPI ${PREMSA} --input ${TMP_FILE} --reference ${REFERENCE} --trim-from ${trim_from} --trim-to ${trim_to} --E 0.01 --N-fraction ${n_frac} --remove-stop-codons Yes
    rm ${TMP_FILE}
    """
}