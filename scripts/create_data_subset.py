#!/usr/bin/env python
import argparse
import pandas as pd
import screed
from pathlib import Path


def create_output_folder(sequences_path, start_date, end_date):
    root_folder = "/".join(sequences_path.split("/")[:-1])
    output_folder = f"{root_folder}/{args.start_date}_{args.end_date}/"
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    return output_folder


def load_metadata(path, date_col):
    df = pd.read_csv(path, low_memory=False).fillna("")
    df.sort_values(by=date_col, inplace=True)
    return df


def shrink_df_by_date(df, start_date, end_date, date_col):
    subset = df[(df[date_col] >= start_date) & (df[date_col] <= end_date)]
    return subset


def separate_random_suspect(df):
    random = df[df["SEQ_REASON"] == "N"]
    suspect = df[df["SEQ_REASON"] != "N"]
    return random, suspect


def load_sequences(path):
    sequences = {}
    with screed.open(path) as seqfile:
        for record in seqfile:
            sequences[record.name] = record.sequence
    return sequences


def write_sequence(output, identifier, sequence):
    line_length = 80
    output.write(f">{identifier}\n")
    formatted_sequence = "".join(
        [
            f"{sequence[i:i + line_length]}\n"
            for i in range(0, len(sequence), line_length)
        ]
    )
    output.write(f"{formatted_sequence}")
    
    
def create_sequence_subset_file(sequences, metadata_subset, output_file, output_folder, seq_id_col):
    with open(output_folder+output_file, "w") as fasta:
        for id in metadata_subset[seq_id_col]:
            if id in sequences:
                write_sequence(fasta, id, sequences[id])
    print("New sequences file:", output_file)
        


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Create a data subset of a specific time period')
    arguments.add_argument('-c', '--csv', required=True, help = 'CSV file containing the metadata', type = str, metavar="PATH")
    arguments.add_argument('-f', '--fasta', required=True, help = 'FASTA file containing the sequence data', type = str, metavar="PATH")
    arguments.add_argument('-s', '--start_date', required=True, help = 'Start date of the time period in the format yyyy-mm-dd', type = str, metavar="DATE")
    arguments.add_argument('-e', '--end_date', required=True, help = 'End date of the time period in the format yyyy-mm-dd', type = str, metavar="DATE")
    arguments.add_argument('-d', '--date_col', required=False, help = 'Name of the column storing the sample date in the metadata [default: DATE_DRAW]', default = "DATE_DRAW", type = str, metavar="STR")
    arguments.add_argument('-i', '--seq_id_col', required=False, help = 'Name of the column storing the sequence identifier in the metadata [default: IMS_ID]', default = "IMS_ID", type = str, metavar="STR")
    arguments.add_argument('--separate', required=True, help = 'Set flag whether or not the output is supposed to be separated into "random" and "suspect" samples', action=argparse.BooleanOptionalAction)
    args = arguments.parse_args()
    
    output_folder = create_output_folder(args.fasta, args.start_date, args.end_date)
    print("Loading metadata...")
    metadata = load_metadata(args.csv, args.date_col)
    print(f"Creating subset based on time frame {args.start_date} to {args.end_date}...")
    subset_time_frame = shrink_df_by_date(metadata, args.start_date, args.end_date, args.date_col)
    del metadata
    if args.separate:
        print("Separating subset into 'random' and 'suspect' samples...")
        subset_random, subset_suspect = separate_random_suspect(subset_time_frame)      
    print("Loading sequences...")
    sequences = load_sequences(args.fasta)
    print(f"Creating output files in folder {output_folder}")
    if args.separate:
        csv_random = f"random.csv"
        csv_suspect = f"suspect.csv"
        subset_random.to_csv(f"{output_folder}{csv_random}", index=False)
        subset_suspect.to_csv(f"{output_folder}{csv_suspect}", index=False)
        print(f"New metadata file: {csv_random}")
        print(f"New metadata file: {csv_suspect}")
        fasta_random = f"random.fasta"
        fasta_suspect = f"suspect.fasta"
        create_sequence_subset_file(sequences, subset_random, fasta_random, output_folder, args.seq_id_col)
        create_sequence_subset_file(sequences, subset_suspect, fasta_suspect, output_folder, args.seq_id_col)
    else:
        csv_subset = f"subset.csv"
        subset_time_frame.to_csv(f"{output_folder}{csv_subset}", index=False)
        print(f"New metadata file: {csv_subset}")
        fasta_subset = f"subset.fasta"
        create_sequence_subset_file(sequences, subset_time_frame, fasta_subset, output_folder, args.seq_id_col)
        