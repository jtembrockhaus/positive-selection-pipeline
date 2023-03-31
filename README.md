# Positive Selection Pipeline

## Requirements
This tool was developed and tested on nextflow (22.04.5) and conda (23.1.0). You can install nextflow using the following conda environment:
```
conda env create -n ps_pipeline -f envs/ps_pipeline.yaml
conda activate ps_pipeline
```
The pipeline required specific HyPhy scripts. Therefore the git repository [HyPhy standalone analyses](https://github.com/veg/hyphy-analyses) need to be cloned to this project and saved in the folder `ressources/hyphy-analyses`.
```
git clone https://github.com/veg/hyphy-analyses.git ressources/hyphy-analyses
```

## Data 
The data to be analyzed has to be in the following format: A file named `sequences.fasta` containing the sequences in Multi-FASTA format and a file named `metadata.csv` containing the sequence information.
The path to the respective directory containing both files nedd to be provided using the `--data_dir` argument.
```
nextflow main.nf --data_dir data/input/desh_subset10
```
Info: If the user wants to create time restricted subsets of the data one can use `scripts/create_data_subset.py` for that purpose.
```
conda env create -n create_data_subset -f envs/create_data_subset.yaml
conda activate create_data_subset
python scripts/create_data_subset.py --help
```