#!/usr/bin/env nextflow

/*
 * 
 */

process SUMMARIZE_SELECTION_ANALYSIS {

    conda "${projectDir}/envs/summarize_selection_analysis.yaml"

    input:
    val gene
    path slac_results_ch
    path fel_results_ch
    path meme_results_ch
    path fubar_results_ch
    path nuc_msa_filtered_ch
    path nuc_msa_variants_duplicates_ch


    output:
    path "${gene}_fubar_results.json", emit: fubar_results_ch
    
    script:
    """
    python ${projectDir}/scripts/summarize_selection_analysis.py \\
    --slac ${slac_results_ch} \\
    --fel ${fel_results_ch} \\
    --meme ${meme_results_ch} \\
    --fubar ${fubar_results_ch}
    --epitopes ${projectDir}/data/static/epitopes.json \\
    %%%--rbd-affinity data/single_mut_effects.csv \\
    %%%--database $MASTERNOFASTA \\
    --coordinates ${nuc_msa_filtered_ch} \\
    --duplicates ${nuc_msa_variants_duplicates_ch} \\
    --pvalue 0.1 \\
    --output summarized_selection_results.json \\
    %%%--evolutionary_annotation data/evo_annotation.json \\
    --mafs ${gene}_mafs.csv \\
    --evolutionary_csv ${gene}_evo_freqs.csv \\
    %%%--evolutionary_fragment $FRAGMENT \\
    %%%--frame_shift ${ADDSHIFT} \\
    %%%--fragment_shift $SHIFT \\
    %%%--offset $OFFSET \\
    --overall $ANNOTATION
    """
}