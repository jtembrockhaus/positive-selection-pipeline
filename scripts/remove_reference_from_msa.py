import sys
import argparse
from Bio import SeqIO

def remove_reference_seq(input_fn, reference_fn, output_fn, removed_seqs_fn):

    seqs = list(SeqIO.parse(input_fn, 'fasta'))
    to_remove = list(SeqIO.parse(reference_fn, 'fasta'))
    to_remove_names = [x.name for x in to_remove]

    to_write = list(filter(lambda x: x.name not in to_remove_names, seqs))
    removed = list(filter(lambda x: x.name in to_remove_names, seqs))

    SeqIO.write(to_write, output_fn, "fasta")
    SeqIO.write(removed, removed_seqs_fn, "fasta")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input',   help = 'fasta to remove sequences from', required = True, type = str)
    arguments.add_argument('-r', '--reference',   help = 'fasta of sequences to remove', required = True, type = str)
    arguments.add_argument('-o', '--output', help = 'fasta without removed sequences', type = str, default = sys.stdout)
    arguments.add_argument('-s', '--removed_seqs', help = 'fasta of removed sequences from msa', type = str, default = sys.stdout)
    args = arguments.parse_args()
    remove_reference_seq(args.input, args.reference, args.output, args.removed_seqs)
