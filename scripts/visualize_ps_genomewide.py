import argparse
import math
import sys
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

gene_lengths = {
    "leader": 540,
    "nsp2": 1914,
    "nsp3": 5835,
    "nsp4": 1500,
    "3C": 918,
    "nsp6": 870,
    "nsp7": 249,
    "nsp8": 594,
    "nsp9":339,
    "nsp10": 417,
    "RdRp":2795,
    "helicase": 1803,
    "exonuclease": 1581,
    "endornase": 1038,
    "methyltransferase": 894,
    "S": 3822,
    "ORF3a": 828,
    "E": 228,
    "M": 669,
    "ORF6": 186,
    "ORF7a": 366,
    "ORF8": 186,
    "N": 908,
    "ORF10": 117
    }

def open_ps_files(paths_str):
    paths = paths_str.split(",")
    
def preprocessing(paths_ps, paths_total):
    sites_under_ps = {}
    total_sites = {}
    for p in paths_ps:
        with open(p, "r") as f:
            gene = p.split("/")[-1].split("_")[0]
            ps = [int(x.replace("\n","")) for x in f.readlines()]
            sites_under_ps[gene] = ps
    for p in paths_total:
        with open(p, "r") as f:
            gene = p.split("/")[-1].split("_")[0]
            total = [int(x.replace("\n","")) for x in f.readlines()]
            total_sites[gene] = total
    print("sites_under_ps:", sites_under_ps)
    print("total_sites:", total_sites)
    
def create_figure():
    fig = px.bar(df, x="Name des Spiels", y="Punktzahl", color="Punktzahl", title="Gesamtbewertung aller Spiele")
    


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='')
    arguments.add_argument('-p', '--sites_ps', required=True, help = 'Comma-separated files with sites under positive selection for a gene', type = str,)
    arguments.add_argument('-t', '--sites_total', required=True, help = 'Comma-separated files with the number of total sites of a gene', type = str,)
    args = arguments.parse_args()
    
    paths_ps = args.sites_ps.split(",")
    paths_total = args.sites_total.split(",")
    
    preprocessing(paths_ps, paths_total)