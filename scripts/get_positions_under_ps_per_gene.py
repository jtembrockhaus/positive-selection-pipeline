import argparse
import math
import sys
import pandas as pd
import numpy as np
from subfuctions import load_json_file
        
def preprocessing_fel(json_dict):    
    df = pd.DataFrame(data=json_dict["MLE"]["content"]["0"])
    df.columns = [x[0] for x in json_dict["MLE"]["headers"]]
    df.insert(loc=0, column="site", value=df.index+1)
    return df

def get_positions_under_PS(df):
    sites_under_ps = []
    for index, row in df.iterrows():
        if row["p-value"] < 0.1:
            sites_under_ps.append(row["site"])
    return sites_under_ps

def write_results(gene, df, sites_under_ps):
    with open(gene + "_fel_sites_under_ps.txt", "w") as f1:
        for site in sites_under_ps:
            f1.write(str(int(site))+"\n")
    with open(gene + "_total_sites.txt", "w") as f2:
        total_positions = df.shape[0]
        f2.write(str(total_positions))
            

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='')
    arguments.add_argument('-g', '--gene', required=True, help = 'Name of the gene', type = str,)
    arguments.add_argument('-f', '--fel', required=True, help = 'FEL analysis results in JSON format', type = argparse.FileType('r'))
    args = arguments.parse_args()
    
    json_dict = load_json_file(args.fel)
    df = preprocessing_fel(json_dict)
    sites_under_ps = get_positions_under_PS(df)
    write_results(args.gene, df, sites_under_ps)
    
    
    