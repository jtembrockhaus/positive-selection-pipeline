import argparse
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from subfuctions import load_json_file, convert_nargs_to_list
from Bio import SeqIO


def get_sample_dates(path, acc_col, date_col):
    df = pd.read_csv(path)
    df = df[[acc_col, date_col]]
    df.rename(columns={acc_col: "accessions", date_col: "date"}, inplace=True)
    df.set_index("accessions", inplace=True)
    return df


def get_file_for_gene(gene, paths_list):
    return [x for x in paths_list if x.split("/")[-1].startswith(gene+"_")][0]


def create_site_subdict(site, gene, position, msa_paths, duplicates_paths):
    sub_dict = {}
    sub_dict["gene"] = gene
    sub_dict["position"] = position
    sub_dict["msa"] = get_file_for_gene(gene, msa_paths)
    sub_dict["duplicates"] = get_file_for_gene(gene, duplicates_paths)
    return sub_dict


def create_sites_dict(genes, sites, msa_paths, duplicates_paths):
    sites_dict = {}
    for site in sites:
        splitted = site.split(":")
        gene = splitted[0]
        position = splitted[1]
        if gene in genes:
            sub_dict = create_site_subdict(site, gene, position, msa_paths, duplicates_paths)
            sites_dict[site] = sub_dict
        else:
            print('No graph could be generated for site "{}" since the gene {} was not analyzed by the pipeline. Change the config to add the gene {} to the pipeline.'.format(site, gene, gene))
    return sites_dict


def add_sites_from_msa(df, msa_path, position):
    df[position] = np.nan
    for record in SeqIO.parse(msa_path, "fasta"):
        acc = record.id.replace("_", "-")
        df.at[acc, position] = record.seq[int(position)-1]
    return df


def check_if_duplicates_exist(duplicates_dict, acc):
    try:
        test = duplicates_dict[acc]["1"]
        return True
    except KeyError:
        return False


def add_sites_for_duplicates(df, duplicates_path, position):
    with open(duplicates_path, "r") as duplicates_file:
        duplicates = load_json_file(duplicates_file)
        for key, sub_dict in duplicates.items():
            has_duplicates = True if len(sub_dict.keys()) > 1 else False
            if has_duplicates:
                for sub_key, value in sub_dict.items():
                    if sub_key != "0":
                        df.at[value.replace("_", "-"), position] = df.at[key.replace("_", "-"), position]
    return df


def get_grouped_df(df, position):
    df = df.groupby(["date", position]).size().reset_index(name='count')
    return df
                
                    
def create_figure(df, gene, position):
    fig = px.bar(df, x="date", y="count", color=position)
    fig.update_yaxes(
        title = "#Sequences with residue",
        titlefont_size = 13,
    )
    fig.update_xaxes(
        title = "Collection date",
        titlefont_size = 13,
    )
    fig.update(
        layout_title = "Temporal Evolution of {}:{}".format(gene, position),
        layout_title_font_size = 24,
        layout_title_x = 0.5,
        layout_legend_title = "Residue"
        )
    fig.write_html("{}_{}_TemporalEvolutionPlot.html".format(gene, position))
        
        

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Create a temporal evolution plot')
    arguments.add_argument('-g', '--genes',  required=True, help = 'List of genes analyzed by the pipeline', type = str, nargs='+')
    arguments.add_argument('-s', '--sites',  required=True, help = 'List of codons to generate plots for', type = str, nargs='+')
    arguments.add_argument('-m', '--msa',  required=True, help = 'Protein Multiple Sequence Alignment of the sequences', type = str, nargs='+')
    arguments.add_argument('-c', '--metadata',  required=True, help = 'Metadata of the sequences', type = argparse.FileType('r'))
    arguments.add_argument('-d', '--duplicates', required=True, help = 'Overview of duplicated sequences', type = str, nargs='+')
    arguments.add_argument('-a', '--acc_col', required=False, help = 'Column name of the metadata file storing the accession numbers', type = str, default="IMS_ID")
    arguments.add_argument('-t', '--date_col', required=False, help = 'Column name of the metadata storing the sampling dates', type = str, default="DATE_DRAW")
    args = arguments.parse_args()
    
    df_metadata = get_sample_dates(args.metadata, args.acc_col, args.date_col)
    
    sites = convert_nargs_to_list(args.sites)
    genes = convert_nargs_to_list(args.genes)
    msa_paths = convert_nargs_to_list(args.msa)
    duplicates_paths = convert_nargs_to_list(args.duplicates)
    
    sites_dict = create_sites_dict(genes, sites, msa_paths, duplicates_paths)
    
    for site, sub_dict in sites_dict.items():
        df_tmp = df_metadata.copy()
        df_site = add_sites_from_msa(df_tmp, sub_dict["msa"], sub_dict["position"])
        df_site = add_sites_for_duplicates(df_site, sub_dict["duplicates"], sub_dict["position"])
        df_site = df_site.dropna() # samples with NaN-values here were filtered out in previous pipeline steps
        df_site = get_grouped_df(df_site, sub_dict["position"])
        create_figure(df_site, sub_dict["gene"], sub_dict["position"])
    