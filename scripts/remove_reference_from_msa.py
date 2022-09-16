import sys
import argparse
from Bio import SeqIO

def remove_reference_seq(input_fn, reference_fn, output_fn):

    seqs = list(SeqIO.parse(input_fn, 'fasta'))
    to_remove = list(SeqIO.parse(reference_fn, 'fasta'))
    to_remove_names = [x.name for x in to_remove]

    to_write = list(filter(lambda x: x.name not in to_remove_names, seqs))

    # Second pass to find nearly similar
    SeqIO.write(to_write, output_fn, "fasta")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input',   help = 'fasta to remove sequences from', required = True, type = str)
    arguments.add_argument('-r', '--reference',   help = 'fasta of sequences to remove', required = True, type = str)
    arguments.add_argument('-o', '--output', help = 'output with removed sequences', type = str, default = sys.stdout)
    args = arguments.parse_args()
    remove_reference_seq(args.input, args.reference, args.output)
