import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import argparse
import warnings
from ete3 import Tree




# read summary file for the NCBI hits
summary_file = f"{snakemake.config['summaries_dir']}/{snakemake.wildcards.cl}.summary"
ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
# remove included=False rows
ncbi_summary_df = ncbi_summary_df[ncbi_summary_df.included]
# replace nan by "NA"
ncbi_summary_df = ncbi_summary_df.fillna("NA")


# read crassvirales taxonomy
crassvirales_taxonomy_file = snakemake.config['crassvirales_taxa_file']
crassvirales_df = pd.read_csv(crassvirales_taxonomy_file, sep="\t", index_col=0, names=["family","subfamily","genus"])


# read the tree
t = Tree(snakemake.input[0], format=1)

# assign taxonomy to the tree
for leaf in t.iter_leaves():
    # check if it is an NCBI hit:
    if leaf.name in ncbi_summary_df.index:
        leaf.add_features(source="ncbi",
                          crass=str(ncbi_summary_df.loc[leaf.name, "ncbi_crass"]),
                          superkingdom=ncbi_summary_df.loc[leaf.name, "superkingdom"],
                          phylum=ncbi_summary_df.loc[leaf.name, "phylum"],
                          #class=ncbi_summary_df.loc[leaf.name, "class"],
                          order=ncbi_summary_df.loc[leaf.name, "order"],
                          family=ncbi_summary_df.loc[leaf.name, "family"],
                          genus=ncbi_summary_df.loc[leaf.name, "genus"]
                          )
    # if not, it is a Crassvirales sequence
    else:
        genome = leaf.name.split("|")[0]
        leaf.add_features(source="crass",
                          crass=str(True),
                          superkingdom="crass",
                          # I put the Crassvirales family here
                          phylum=crassvirales_taxa_df.loc[genome, "family"],
                          #class="crass",
                          order="crass",
                          family=crassvirales_taxa_df.loc[genome, "family"],
                          genus=crassvirales_taxa_df.loc[genome, "genus"]
                          )


# identify all the monophyletic clades containing Crassvirales, Viruses or NA
crassvirales_clades = t.get_monophyletic(values=["crass", "new", "Viruses", "NA"], target_attr="superkingdom")
