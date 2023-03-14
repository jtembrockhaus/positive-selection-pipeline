import argparse
import os
from subfuctions import load_json_file
from Bio import Phylo
from Bio.Phylo.PAML import codeml

def analyze_duplicates(duplicates_json):
    number_of_duplicates = 0
    max_identical_seqs = 0
    for seq in duplicates_json:
        dups = len(duplicates_json[seq].keys())
        number_of_duplicates += dups if dups > 1 else 0
        max_identical_seqs = dups if dups > max_identical_seqs else max_identical_seqs
    return number_of_duplicates, max_identical_seqs

def convert_to_relative_path(path, root_dir):
    return path.replace(root_dir+"/","")

# def analyze_tree(gene, pipeline_dir, tree, msa, jones):
#     cml = codeml.Codeml(
#     alignment=convert_to_relative_path(msa, pipeline_dir),
#     tree=convert_to_relative_path(tree, pipeline_dir),
#     out_file=convert_to_relative_path(os.getcwd()+"/"+gene+"_tree_analysis.json", pipeline_dir),
#     working_dir=pipeline_dir,
#     )
#     cml.set_options(aaRatefile="/".join([pipeline_dir, jones]))
#     cml.run(verbose=True)

#TODO Give statistics for phylogenetic tree
#TODO Give statistics for filtered samples (low frequency variants, n-content, ...)

def create_report(number_of_duplicates, max_identical_seqs, output_path):
    with open(output_path, "w") as out_file:
        out_file.write("Number of duplicates: " + str(number_of_duplicates) + "\n")
        out_file.write("Max. Number of identical sequences: " + str(max_identical_seqs) + "\n")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Create a report file for the pipeline run of a single gene')
    arguments.add_argument('-g', '--gene', help = 'Name of the gene', type = str,)
    arguments.add_argument('-p', '--pipeline_dir',   help = 'Root Directory of the pipeline', type = str)
    arguments.add_argument('-d', '--duplicated_seqs',   help = 'Overview of duplicated sequences in JSON format', type = argparse.FileType('r'))
    arguments.add_argument('-t', '--tree',   help = 'Phylogenetic tree in Newick format', type = str)
    arguments.add_argument('-a', '--msa',   help = 'Multiple Sequence Alignment (MSA) on nucleotide level in FASTA format', type = str)
    arguments.add_argument('-j', '--jones',   help = 'jones.dat file of the Codeml package [default: data/static/codeml/dat/jones.dat]', default = "data/static/codeml/dat/jones.dat", type = str)
    arguments.add_argument('-o', '--output_file', help = 'Path of the output report file', type = str,)
    args = arguments.parse_args()
    
    duplicates_json = load_json_file(args.duplicated_seqs)
    number_of_duplicates, max_identical_seqs = analyze_duplicates(duplicates_json)
    #analyze_tree(args.gene, args.pipeline_dir, args.tree, args.msa, args.jones)
    create_report(number_of_duplicates, max_identical_seqs, args.output_file)
