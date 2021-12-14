import pandas as pd
import numpy as np

def set_vOTUs(votus_file):
    '''
    Reads the two-columns tabular file with the vOTU representative and the rest of genomes
    in the vOTU cluster. Store in a dict, k=genome v=representative.

    Return a second dictionary with key=representative  v={size:0, votu:[genome1, genome2, genome3...]}
    '''

    lines = [line.strip().split("\t") for line in open(votus_file).readlines()]

    reprs_votus = {line[0]:line[1].split(",") for line in lines}

    genomes_reprs = dict()
    for line in lines:
        genomes = line[1].split(",")
        for genome in genomes:
            genomes_reprs[genome] = line[0]


    votus = {line[0]:{"size":len(line[1].split(",")), "votu":line[1].split(",")} for line in lines}

    return genomes_reprs, reprs_votus, votus

def get_taxonomy_df(taxa_file):
    '''
    Parses the taxonomy and returns a dataframe, index=genomes, 3 columns (family, subfamily, genus)
    '''
    #"/home/danielc/software/github/old_june/tree_functions/files/crassphages_taxonomy_terL_and_new.txt"
    df = pd.read_csv(taxa_file, sep="\t", index_col=0, names=["family","subfamily","genus"])
    df['superkingdom'] = np.where(df['family'].isin(["new", "crass_env"]), "new", "crassvirales")
    return df
