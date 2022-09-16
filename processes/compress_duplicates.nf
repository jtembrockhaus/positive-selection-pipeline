#!/usr/bin/env nextflow

/*
 * Create unique sets of protein and nucleotide haplotypes  
 * while keeping track of what other sequences are being represented
 * by the retained unique sequences.
 */

process COMPRESS_DUPLICATES {

    conda "${projectDir}/envs/compress_duplicates.yaml"

    input:
    val gene
    path gene_protein_ch
    path gene_nuc_ch

    output:
    path "${gene}_protein_compressed.fas", emit: gene_protein_compressed_ch
    path "${gene}_nuc_compressed.fas", emit: gene_nuc_compressed_ch
    path "${gene}_protein_duplicates.json", emit: gene_protein_duplicates_ch
    path "${gene}_nuc_duplicates.json", emit: gene_nuc_duplicates_ch

    script:
    """
    python ${projectDir}/scripts/compress_duplicates.py \
    --protein-input ${gene_protein_ch} \
    --nuc-input ${gene_nuc_ch} \
    --protein-output ${gene}_protein_compressed.fas \
    --nuc-output ${gene}_nuc_compressed.fas \
    --protein-duplicates ${gene}_protein_duplicates.json \
    --nuc-duplicates ${gene}_nuc_duplicates.json
    """
}