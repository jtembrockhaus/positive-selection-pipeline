import csv
import json
import sys
import argparse
import copy
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Create unique sets of protein and nucleotide haplotypes while keeping track of what other sequences are being represented by the retained unique sequences')
arguments.add_argument('--protein-input', help = 'protein fasta to filter for duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('--nuc-input', help='nucleotide fasta to filter for duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('--protein-output', help = 'compressed protein fasta output file', required = True, type = argparse.FileType('w'))
arguments.add_argument('--nuc-output', help = 'compressed nucleotide fasta output file', required = True, type = argparse.FileType('w'))
arguments.add_argument('--protein-duplicates', help='protein duplicates lookup JSON output file', required = True, type = argparse.FileType('w'))
arguments.add_argument('--nuc-duplicates', help='nucleotide protein duplicates lookup JSON output file', required = True, type = argparse.FileType('w'))
args = arguments.parse_args()


def get_duplicates(seqs):
    seqs = sorted(seqs, key=lambda x: x.seq)
    duplicates = {}
    for rec in seqs:
        if(rec.seq not in duplicates):  
            duplicates[rec.seq] = [rec]
        else:
            duplicates[rec.seq].append(rec)
    return duplicates


protein_seqs = list(SeqIO.parse(args.protein_input, 'fasta'))
protein_duplicates = get_duplicates(protein_seqs)

nuc_seqs = list(SeqIO.parse(args.nuc_input, 'fasta'))
nuc_duplicates = get_duplicates(nuc_seqs)

# Write new fasta from first index of each dupe, and duplicates file if provided
vals = list(protein_duplicates.values())

# Write fasta
firsts = [val[0] for val in vals]

# Write protein_duplicates
dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in vals}
json.dump(dupe_names, args.protein_duplicates, indent=4, sort_keys=True)

# Write nuc_duplicates
nuc_vals = list(nuc_duplicates.values())
nuc_dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in nuc_vals}
json.dump(nuc_dupe_names, args.nuc_duplicates, indent=4, sort_keys=True)

for first in firsts:
    first.id = first.name
    first.description = first.name
    # first.id = first.name + "_" + str(len(dupe_names[first.name].keys()))
    # first.description = first.name + "_" + str(len(dupe_names[first.name].keys()))

nuc_firsts = [val[0] for val in nuc_vals]
filtered_seq_names = [seq.name for seq in nuc_firsts]

for first in nuc_firsts:
    first.id = first.name
    first.description = first.name
    # first.id = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))
    # first.description = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))


# Protein seq map
prot_seq_map = {record.id : record for record in protein_seqs}
prot_nuc_firsts = [prot_seq_map[nuc_first.id] for nuc_first in nuc_firsts]

# Write protein sequences based on nucleotide dupes
# Get protein_seqs based on nucleotide dupes
SeqIO.write(prot_nuc_firsts, args.protein_output, "fasta")
SeqIO.write(nuc_firsts, args.nuc_output, "fasta")