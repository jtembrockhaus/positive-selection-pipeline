#!/usr/bin/env nextflow

/*
 * Estimate the probability that singleton low frequency variants
 * are errors and replace them with "-". Low frequency variants
 * that co-occur are not removed. Sequences that have too many low
 * frequency variants are removed.
 */

process FILTER_NUC_MSA {

    conda "${projectDir}/envs/filter_nuc_msa.yaml"

    input:
    val gene
    path nuc_msa_merged_ch
    path nuc_msa_merged_duplicates_ch

    output:
    path "${gene}_nuc_msa_filtered.fas", emit: nuc_msa_filtered_ch
    

    script:
    COMPRESSOR="compressor.bf"
    COMPRESSOR2="compressor_2.bf"
    """
    cp ${projectDir}/scripts/${COMPRESSOR} .
    cp ${projectDir}/scripts/${COMPRESSOR2} .
    hyphy ${COMPRESSOR} --msa ${nuc_msa_merged_ch} --duplicates ${nuc_msa_merged_duplicates_ch} --output ${gene}_nuc_msa_variants.csv --json ${gene}_nuc_msa_variants.json --duplicate-out ${gene}_nuc_msa_variants_duplicates.json
    hyphy ${COMPRESSOR2} --msa ${nuc_msa_merged_ch} --duplicates ${nuc_msa_merged_duplicates_ch} --csv ${gene}_nuc_msa_variants.csv --byseq ${gene}_nuc_msa_variants.json --p 0.95 --output ${gene}_nuc_msa_filtered.fas --json ${gene}_nuc_msa_filtered.json --output-edits ${gene}_nuc_msa_filtered_edits.json
    rm ${COMPRESSOR}
    rm ${COMPRESSOR2}
    """
}