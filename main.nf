#!/usr/bin/env nextflow

/*
Nextflow -- Analysis pipeline to determine positive selection in SARS-CoV-2 samples
Author: Julius Tembrockhaus
Based on a research project by Sergei Pond
*/


/**************************
* HELP MESSAGE
**************************/
if (params.help) { exit 0, help_message() }
// TODO Add a check whether hyphy-analysis was already set up


/**************************
* INITIALIZATIONS
**************************/
import groovy.json.JsonSlurper


/**************************
* INPUT ARGUMENTS
**************************/
// --data_dir
if (!params.data_dir){
  exit 1, "No data directory specified. Please provide a folder containing a file named sequences.fasta and metadata.csv using the --data_dir argument!"
}
data_dir = params.data_dir.replaceAll("/\\z", "") // remove trailing slash if it exist
sequences_file = new File("$data_dir/sequences.fasta")
metadata_file = new File("$data_dir/metadata.csv")
if (!sequences_file.exists() || !metadata_file.exists()) {
  exit 1, "The provided data directory $params.data_dir must contain a file called sequences.fasta and a file called metadata.csv!"
} else {
  sequences_ch = Channel.value("$projectDir/$data_dir/sequences.fasta")
  metadata_ch = Channel.value("$projectDir/$data_dir/metadata.csv")
}
// TODO Add a check for required columns in metadata
// TODO Add the possibility to use compressed input files (e.g. sequences.fasta.gz)
// --gene_list
gene_file = params.gene_list ? new File(params.gene_list) : new File("configs/genes_and_products.txt")
if (!gene_file.exists()){
  exit 1, "File $params.gene_list provided for --gene_list does not exist."
} else {
  genes = get_genes_from_file(gene_file)
  genes_ch = Channel.from(genes)
}
// --gene_trim_interval
trim_file = params.gene_trim_interval ? new File(params.gene_trim_interval) : new File("data/static/padded_gene_intervals.json")
if (!trim_file.exists()){
  exit 1, "File $params.gene_trim_interval provided for --gene_trim_interval does not exist."
} else {
  trim_intervals = get_gene_trim_intervals(trim_file, genes)
  trim_from_ch = Channel.from(trim_intervals*.get(0))
  trim_to_ch = Channel.from(trim_intervals*.get(1))
  n_frac_ch = Channel.from(trim_intervals*.get(2))
}

/**************************
* PROCESSES
**************************/
include { EXTRACT_GENE } from "./processes/extract_gene.nf"


/**************************
* MAIN WORKFLOW 
**************************/
workflow {
  EXTRACT_GENE(sequences_ch, metadata_ch, genes_ch, trim_from_ch, trim_to_ch, n_frac_ch)
  gene_protein_ch = EXTRACT_GENE.out.gene_protein_ch
  gene_nuc_ch = EXTRACT_GENE.out.gene_nuc_ch
}

/*************
* --help
*************/
def help_message() {
    log.info """
    ____________________________________________________________________________________________
    Workflow: Positive Selection Analysis Pipeline

    Usage example:
    nextflow run main.nf 

    Mandatory arguments:
    --data_dir                  Path to directory containing the sequences in FASTA format and the metadata in CSV format

    Optional arguments:
    --gene_list                 Path to file containing a list of genes and products to be analyzed
                                [default: configs/genes_and_products.txt]
    --gene_trim_interval        Path to file containing padded coordinate-ranges where the gene is expected to be 
                                and tolerated N fractions for the gene extraction in JSON format
                                [default: data/static/padded_gene_intervals.json]

    Note: Paths of listed folders need to be relative to location of main.nf

    """.stripIndent()
}


/**************************
* FUNCTIONS
**************************/
def get_genes_from_file(file){
  genes = file.readLines()
  genes.removeAll{it.startsWith('#')}
  return genes
}
def get_gene_trim_intervals(file, genes){
  def jsonSlurper = new JsonSlurper()
  data = jsonSlurper.parse(file)
  def gene_trim_intervals = []
  for (gene in genes) {
    gene_trim_intervals << new Tuple(
      data[gene]["trim_from"],
      data[gene]["trim_to"],
      data[gene]["n_frac"]
      )
  }
  return gene_trim_intervals
}