#!/usr/bin/env nextflow

/*
 * Merge duplicates from post-msa and prior.
 */

process MERGE_DUPLICATES {

    conda "${projectDir}/envs/merge_duplicates.yaml"

    input:
    val gene
    path nuc_duplicates_ch
    path nuc_msa_duplicates_ch
    path nuc_msa_compressed_ch

    output:
    path "${gene}_nuc_msa_merged_duplicates.fas", emit: nuc_msa_merged_duplicates_ch

    script:
    MERGED_DUPLICATES="${gene}_nuc_merged_duplicates.json"
    MAPPING="${gene}_nuc_merged_duplicates_map.json"
    """
    python ${projectDir}/scripts/merge_duplicates.py -p ${nuc_duplicates_ch} -n ${nuc_msa_duplicates_ch} -o ${MERGED_DUPLICATES}
    python ${projectDir}/scripts/fix_duplicates.py -d ${MERGED_DUPLICATES} -m ${MAPPING} -o
    python ${projectDir}/scripts/update_fasta_duplicates.py -f ${nuc_msa_compressed_ch} -m ${MAPPING} -o ${gene}_nuc_msa_merged_duplicates.fas
    """
}