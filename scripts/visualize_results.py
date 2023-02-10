import argparse
import json
import math
import sys
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff


# CSS color names for classes
classes_color_dict = {
    "Purifying": "palegreen",
    "Diversifying": "palevioletred",
    "Neutral": "lightslategray",
    "Invariable": "white"
}

# dictionary for titles of MEME variable plots
meme_variable_dict = {
    'α': "Synonymous substitution rate",
    'β⁻': "Non-synonymous substitution rate for the negative/neutral evolution component",
    'p⁻': "Mixture distribution weight allocated to β⁻",
    'β⁺': "Non-synonymous substitution rate for the positive/neutral evolution component",
    'p⁺': "Mixture distribution weight allocated to β⁺",
    'LRT': "Likelihood ratio test statistic for episodic diversification",
    'p-value': "Asymptotic p-value for episodic diversification",
    '# branches under selection': "Approximate estimate of how many branches may have been under selection",
    'Total branch length': "Total length of branches contributing to inference and used to scale dN-dS",    
}

# function to convert to superscript
def get_super(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = x.maketrans(''.join(normal), ''.join(super_s))
    return x.translate(res)
    
# load FEL results
def load_json_file(filename):
    try:
        json_dict = json.load(filename)
        return json_dict
    except Exception as e:
        print(e)
        sys.exit(1)

# censor values at 10 such that values >10 are replaced by 10      
def censoring_at_10(value):
    if value == np.nan:
        return np.nan
    elif value <= 10:
        return value
    else:
        return 10
        
# preprocess FEL data
def preprocessing_fel(json_dict):
    df = pd.DataFrame(data=json_dict["MLE"]["content"]["0"])
    df.columns = [x[0] for x in json_dict["MLE"]["headers"]]
    df.insert(loc=0, column="site", value=df.index+1)
    df["class"] = "unassigned"
    for index, row in df.iterrows():
        if row["p-value"] == 1:
            df.at[index, "class"] = "Invariable"
        elif row["p-value"] > 0.1:
            df.at[index, "class"] = "Neutral"
        elif row["alpha"] > row["beta"]:
            df.at[index, "class"] = "Purifying"
        elif row["beta"] > row["alpha"]:
            df.at[index, "class"] = "Diversifying"
        else:
            print("Class cannot be determined for:\n", row) 
    df["alpha=beta_base"] = df["alpha=beta"].apply(lambda x: x/-2)
    df["alpha_censored"] = df["alpha"].apply(lambda x: censoring_at_10(x))
    df["beta_censored"] = df["beta"].apply(lambda x: censoring_at_10(x))
    df["dNdS"] = df["beta"] / df["alpha"]
    df["dNdS"].replace([np.inf, -np.inf], np.nan, inplace=True) # inf occurs when diving x by 0
    df["dNdS_censored"] = df["dNdS"].apply(lambda x: censoring_at_10(x))
    df["alpha"] = df["alpha"].apply(lambda x: x*-1)
    return df

# preprocess MEME data
def preprocessing_meme(json_dict):
    df = pd.DataFrame(data=json_dict["MLE"]["content"]["0"])
    df.columns = [x[0] for x in json_dict["MLE"]["headers"]]
    df.insert(loc=0, column="site", value=df.index+1)
    df.rename(columns = {
        '&alpha;': 'α',
        '&beta;<sup>-</sup>': 'β{}'.format(get_super('-')),
        'p<sup>-</sup>': 'p{}'.format(get_super('-')),
        '&beta;<sup>+</sup>': 'β{}'.format(get_super('+')),
        'p<sup>+</sup>': 'p{}'.format(get_super('+')),        
        }, inplace = True)
    return df

# initialize desired subplots for MLE figure
def setup_subplots_mle(df, codons_per_subplot):
    codons = df.shape[0]
    number_subplots = math.ceil(codons/codons_per_subplot)
    fig = make_subplots(rows=number_subplots, cols=1)
    return fig, codons, number_subplots

# create individual subset
def create_subset(df, i, codons, number_subplots, codons_per_subplot):
    if i == number_subplots:
        start = (i-1)*codons_per_subplot+1
        end = codons-1
    else:
        start = (i-1)*codons_per_subplot
        end = i*codons_per_subplot
    df_subset = df[start:end] 
    return df_subset

# assign colors to class values
def assign_colors(key):
    return classes_color_dict[key]
    
# manually setup legend for each class in the MLE figure
def setup_mle_legend(fig, classes_color_dict):
    for key in classes_color_dict:
        fig.add_trace(go.Bar(
            x=[None], y=[None],
            name=key,
            marker=dict(color=classes_color_dict[key])
            )
        )
    return fig
    
# create figure for MLE visualization of FEL results
def create_mle_figure(df, fig, codons, number_subplots, codons_per_subplot):
    for i in range(1, number_subplots+1):
        df_subset = create_subset(df, i, codons, number_subplots, codons_per_subplot)
        # MLE synonymous rate (α) and non-synonymous rate (β)
        for col in ["alpha", "beta"]:
            fig.add_trace(
                go.Bar(
                    x=df_subset["site"],
                    y=df_subset[col],
                    marker_color=list(map(assign_colors, df_subset["class"])),
                    showlegend=False,
                ),
                row=i, col=1,
            )
        # Estimates under the null model (α=β)
        fig.add_trace(
            go.Bar(
                x=df_subset["site"],
                y=df_subset["alpha=beta"],
                base=df_subset["alpha=beta_base"],
                marker_color='dimgray',
                name='alpha=beta',
                width=0.2,
                showlegend=False,
            ),
            row=i, col=1,
        )
    # initialize legend
    fig = setup_mle_legend(fig, classes_color_dict)
    # x-axis customizations
    fig.update_xaxes(
        title = "Codon",
        titlefont_size = 13,
        dtick = 1,
        tickfont_size = 9,
        tickangle = 270,
    )
    # y-axis customization
    fig.update_yaxes(
        title = "Rate estimate",
        titlefont_size = 13,
        range = [-10,10],
        tickvals = [-10, -5, 0, 5, 10],
        ticktext = ["α=10", "α=5", "0", "β=5", "β=10"],
        tickfont_size = 9,
    )
    # general layout customizations
    fig.update(
        layout_barmode = "overlay",
        layout_legend_title = "Classes",
        layout_legend = dict(font = dict(size = 13)),
        layout_title = "MLE of synonymous (α) and non-synonymous (β) rates at each site",
        layout_title_font_size = 24,
        layout_title_x = 0.5,
        layout_autosize = True,
        layout_height = 250*number_subplots,
        )
    return fig

# add subplot for specific substitution rate into rate density figure
def add_rd_subplot(fig, df, col, name, subplot_position):
    tmp_fig = ff.create_distplot([df[col].values], [name])
    fig.add_trace(
        go.Scatter(
            tmp_fig["data"][1],
            line=dict(color='dimgray'),
            fill='tozeroy',
            ),
        row=subplot_position, col=1,
    )
    mean = abs(np.mean(df[col.replace("_censored","")])) # unconsored 
    #mean2 = abs(np.mean(df[col])) #censored
    fig.add_shape(go.layout.Shape(type="line", x0=mean, x1=mean, y0=0, y1=1, line=dict(color="firebrick")), row=subplot_position, col=1)
    fig.add_annotation(x=mean+0.1, y=0.92, text=str(round(mean,2)), font=dict(size=9, color="gray"), showarrow=False, row=subplot_position, col=1)
    fig.update_xaxes(
        title = name,
        titlefont_size = 13,
        dtick = 0.5,
        tickfont_size = 9,
        row = subplot_position,
        col = 1,
    )
    fig.update_yaxes(
        title = "Density",
        titlefont_size = 13,
        tickfont_size = 9,
        dtick = 0.2,
    )
    return fig

# create figure for rate density visualization of FEL results
def create_rd_figure(df):
    fig = make_subplots(rows=3, cols=1)
    fig = add_rd_subplot(fig, df, "alpha_censored", "α", 1)
    fig = add_rd_subplot(fig, df, "beta_censored", "β", 2)
    fig = add_rd_subplot(fig, df, "dNdS_censored", "dN/dS", 3)
    fig.update(
        layout_title = "Kernel density estimates of site-level rate estimates",
        layout_title_font_size = 24,
        layout_title_x = 0.5,
        layout_autosize = True,
        layout_height = 300*3,
        )
    return fig 

# create histogram figure for one specific variable
def create_histogram_figure(df, fig, variable, codons, number_subplots, codons_per_subplot):
    for i in range(1, number_subplots+1):
        df_subset = create_subset(df, i, codons, number_subplots, codons_per_subplot)
        fig.add_trace(
            go.Bar(
                x=df_subset["site"],
                y=df_subset[variable],
                marker_color="cornflowerblue",
                showlegend=False,
            ),
            row=i, col=1,
        )
    # x-axis customizations
    fig.update_xaxes(
        title = "Codon",
        titlefont_size = 13,
        dtick = 50,
        tickfont_size = 9,
    )
    # y-axis customization
    fig.update_yaxes(
        title = variable,
        titlefont_size = 13,
        tickfont_size = 9,
    )
    # general layout customizations
    fig.update(
        layout_title = meme_variable_dict[variable]+ " at each site",
        layout_title_font_size = 24,
        layout_title_x = 0.5,
        layout_autosize = True,
        layout_height = 250*number_subplots,
        )
    return fig 

# save data frame with desired columns as TSV file   
def save_fel_df(df, gene):
    columns_to_keep = ["site", "alpha", "beta", "alpha=beta", "LRT",
                       "p-value", "Total branch length", "p-asmp", "class"]
    df = df[columns_to_keep]
    df["alpha"] = df["alpha"].apply(lambda x: x*-1)
    df.to_csv(gene+"_FEL.tsv", sep="\t", index=False)

# create MLE plot and save it into HTML file            
def save_mle_figure(df, codons_per_subplot, gene):
    fig, codons, number_subplots = setup_subplots_mle(df, codons_per_subplot)
    fig = create_mle_figure(df, fig, codons, number_subplots, codons_per_subplot)
    fig.write_html(gene+"_FEL_mle.html")
    
# create rate density plot and save it into HTML file    
def save_rd_figure(df, gene):
    fig = create_rd_figure(df)
    fig.write_html(gene+"_FEL_rate_densities.html")
    
# create a histogram for each desired variable and sve it into HTML file
def save_meme_figures(df, codons_per_subplot, gene):
    desired_variables = ['α', 'β⁻', 'p⁻', 'β⁺', 'p⁺', 'LRT', 'p-value','# branches under selection', 'Total branch length']
    for var in desired_variables:
        fig, codons, number_subplots = setup_subplots_mle(df, codons_per_subplot)
        fig = create_histogram_figure(df, fig, var, codons, number_subplots, codons_per_subplot)
        fig.write_html(gene+"_MEME_" + var + ".html")
     

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Create plot(s) based on selection analysis results')
    arguments.add_argument('-g', '--gene',   help = 'Name of the gene', type = str,)
    arguments.add_argument('-f', '--fel',   help = 'FEL analysis results in JSON format', type = argparse.FileType('r'))
    arguments.add_argument('-m', '--meme',   help = 'MEME analysis results in JSON format', type = argparse.FileType('r'))
    args = arguments.parse_args()
    
    if args.fel:
        json_dict = load_json_file(args.fel)
        df = preprocessing_fel(json_dict)
        codons_per_subplot = 100
        save_fel_df(df, args.gene)
        save_mle_figure(df, codons_per_subplot, args.gene)
        save_rd_figure(df, args.gene)
        
    if args.meme:
        json_dict = load_json_file(args.meme)
        df = preprocessing_meme(json_dict)
        codons_per_subplot = 500
        save_meme_figures(df, codons_per_subplot, args.gene)