import argparse
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from subfuctions import load_json_file

    
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
            total = int(f.readline())
            total_sites[gene] = total
    return sites_under_ps, total_sites
    
def create_figure(sites_under_ps, total_sites, gene_lengths):
    fig = px.bar()
    start_pos = 0
    box_height = 5.5
    elevation = 0.5
    gene_starts = {}
    elevated = {}
    # setup every 2nd gene to be elevated for better visualization
    for gene_name in gene_lengths.keys():
        elevated[gene_name] = True if list(gene_lengths.keys()).index(gene_name) % 2 else False
    # add the shape for each gene to the plot
    for gene_name, ref_gene_length in gene_lengths.items():
        num_codons = total_sites[gene_name] if gene_name in total_sites else ref_gene_length
        fig.add_shape(type="rect",
                    x0=start_pos,
                    x1=start_pos+num_codons,
                    y0=elevation if elevated[gene_name] else 0,
                    y1=box_height+elevation if elevated[gene_name] else box_height,
                    line=dict(color="grey",
                              width=2),
                    opacity=1 if gene_name in sites_under_ps.keys() else 0.4,
                    fillcolor="white" if gene_name in sites_under_ps.keys() else "ghostwhite",
                    )
        fig.add_trace(go.Scatter(
            x=[start_pos+num_codons*0.5],
            y=[box_height+elevation+0.5] if elevated[gene_name] else [-0.5],
            text=[gene_name],
            mode="text",
        ))
        gene_starts[gene_name] = start_pos
        start_pos = start_pos + num_codons
    # mark the positive selected sites for every analyzed gene
    for gene_name, sites in sites_under_ps.items():
        for s in sites:
            fig.add_shape(type="line",
                          x0=gene_starts[gene_name]+s,
                          x1=gene_starts[gene_name]+s,
                          y0=0.7 if elevated[gene_name] else 0.2,
                          y1=box_height+elevation-0.2 if elevated[gene_name] else box_height-0.2,
                          line=dict(color="red",
                                    width=0.5)
                          )
    # add a custom legend entry
    fig.add_shape(type="rect",
                  x0=start_pos*0.8,
                  x1=start_pos*0.81,
                  y0=9*0.99,
                  y1=9*1.015,
                  fillcolor="red",
                  line=dict(color="red"))
    fig.add_trace(go.Scatter(
        x=[start_pos-1000],
        y=[9],
        text=[" = Positively selected"],
        mode="text",
        textposition="middle left",
        textfont=dict(size=15)
    ))
    fig.update_yaxes(
        range=[-2,10],
        showticklabels=False,
    )
    fig.update_xaxes(
        title = "Codons",
        titlefont_size = 13,
    )
    fig.update(
        layout_title = "Genomewide Positive Selection",
        layout_title_font_size = 24,
        layout_title_x = 0.5,
        layout_font=dict(size=10), # size of gene labels
        layout_showlegend = False,
        )
    # fig.show()
    return fig

def save_figure(fig):
    fig.write_html("fel_genomewide_positive_selection.html")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Create plot to visualize positively selected sites in the entire genome')
    arguments.add_argument('-p', '--sites_ps', required=True, help = 'One or more files with sites under positive selection for a gene', type = str, nargs='+')
    arguments.add_argument('-t', '--sites_total', required=True, help = 'One or more files with the number of total sites of a gene', type = str, nargs='+')
    arguments.add_argument('-l', '--gene_lengths', required=True, help = 'JSON file containing gene lengths', type = argparse.FileType('r'))
    args = arguments.parse_args()
    
    # paths_ps = args.sites_ps.replace("[","").replace(", ",",").split(",")
    # paths_total = args.sites_total.replace("[","").replace(", ",",").split(",")
    
    paths_ps = [x.replace("[","").replace("]","").replace(",","") for x in args.sites_ps]
    paths_total = [x.replace("[","").replace("]","").replace(",","") for x in args.sites_total]
    
    gene_lengths = load_json_file(args.gene_lengths)
    sites_under_ps, total_sites = preprocessing(paths_ps, paths_total)
    fig = create_figure(sites_under_ps, total_sites, gene_lengths)
    save_figure(fig)