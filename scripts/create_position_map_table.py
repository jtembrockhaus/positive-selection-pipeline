from Bio import SeqIO
import pandas as pd
import argparse


def create_position_map_table(ref_from_msa_fn, pos_map_table_fn):
    mapped_reference = SeqIO.read(ref_from_msa_fn, "fasta").seq
    pos = 1
    suffix = 1
    pos_table = pd.DataFrame(columns=["reference", "msa"])
    for i, base in enumerate(mapped_reference):
        if base != "-": 
            row = pd.DataFrame([[str(pos), i+1]],columns=["reference", "msa"])
            pos_table = pd.concat([pos_table, pd.DataFrame(row, columns=["reference", "msa"])])
            suffix = 1
            pos += 1

        else:
            row = pd.DataFrame([[str(pos-1) + "." + str(suffix), i+1]],columns=["reference", "msa"])
            pos_table = pd.concat([pos_table, pd.DataFrame(row, columns=["reference", "msa"])])
            suffix += 1

    pos_table.to_csv(pos_map_table_fn, sep="\t", index=False)

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='creates position table for convinient numbering of bases')
    arguments.add_argument('-r', '--ref_from_msa', type=str, help='Input FASTA file containing reference sequence mapped to MSA')
    arguments.add_argument('-p', '--pos_map_table', type=str, help='Output TSV file containing positions of reference and MSA')
    args = arguments.parse_args()
    
    create_position_map_table(args.ref_from_msa, args.pos_map_table)

