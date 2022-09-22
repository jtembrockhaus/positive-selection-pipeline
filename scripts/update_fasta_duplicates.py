import json
import argparse
from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-f', '--fasta-file',   help = 'fasta to overwrite', required = True, type = argparse.FileType('r'))
arguments.add_argument('-m', '--map-file',   help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'write updated fasta here', required = True, type = argparse.FileType('w'))
args = arguments.parse_args()

# If one fails, then copy the other to the output. If both fail, then throw an error
map_json = json.load(args.map_file)
seqs = list(SeqIO.parse(args.fasta_file, 'fasta'))

# Fix FASTA headers
for seq in seqs:
    old_id = seq.id
    seq.id = map_json[old_id]
    seq.description = map_json[old_id]

args.fasta_file.close()

SeqIO.write(seqs, args.output, "fasta")
args.output.close()