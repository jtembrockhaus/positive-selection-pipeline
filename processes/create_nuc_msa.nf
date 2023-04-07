#!/usr/bin/env nextflow

/*
 * Map frameshift corrected nucleotide codon sequences onto 
 * protein MSA to obtain a nucleotide MSA. The output is 
 * also compressed to unique sequences while keeping track of
 * other sequences are being represented by the retained 
 * unique sequences. Finally, unknown characters in the output
 * are replaced with N.
 */

process CREATE_NUC_MSA {

    conda "${projectDir}/envs/extract_gene.yaml"

    input:
    val gene
    path protein_msa_ch
    path nuc_seqs_ch

    output:
    path "${gene}_nuc_msa_compressed.fas", emit: nuc_msa_compressed_ch
    path "${gene}_nuc_msa_duplicates.json", emit: nuc_msa_duplicates_ch

    script:
    POSTMSA="${projectDir}/ressources/adjusted-hyphy-analyses/post-msa.bf"
    """
    hyphy ${POSTMSA} --protein-msa ${protein_msa_ch} --nucleotide-sequences ${nuc_seqs_ch} --output ${gene}_nuc_msa_compressed.fas --duplicates ${gene}_nuc_msa_duplicates.json --compress Yes
    sed -i '/^>/! s/[^ACTG-]/N/g' ${gene}_nuc_msa_compressed.fas
    """
}