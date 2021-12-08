'''
It adds full rank information to the summary table of the NR hits.
'''

from ete3 import NCBITaxa
import os, glob, argparse, multiprocessing, time
from functools import partial
import pandas as pd
import numpy as np



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--summaries_dir',
                               dest='summaries_dir',
                               required=True,
                               help='directory with the summary tables obtained with filter_nr_hits.py. '
                               'Output tables will be in the same dir with a "_ranks" in the extension.'
                               )

    return parser.parse_args()

def timeit(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

def add_lineage_to_table(summary_table):
    ncbi = NCBITaxa()
    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus"]

    columns = ["protein_id", "included", "evalue", "bitscore", "p_len", "cl_len", "p_cov", "cl_cov", "lca", "lca_rank", "n_genomes"]
    table = pd.read_csv(summary_table, names=columns, index_col=0, sep="\t", low_memory=False)

    # count included sequences
    n_included = table.included.sum()
    if n_included > 0:
        # filter table df to keep only included
        included = table.loc[table.included, :]

        for hit in included.index:
            # init the ranks as "NA"
            table.loc[hit, "superkingdom"] = "not_resolved"
            table.loc[hit, "phylum"] = "not_resolved"
            table.loc[hit, "order"] = "not_resolved"
            table.loc[hit, "family"] = "not_resolved"
            table.loc[hit, "genus"] = "not_resolved"

            # proces
            taxa_name = included.loc[hit, "lca"]
            try:
                taxid = ncbi.get_name_translator([taxa_name])[taxa_name][0]
                lineage = ncbi.get_lineage(taxid)
                ranks = ncbi.get_rank(lineage)

            except:







def main():

    args = parse_args()

    # list summary tables
    summary_tables = glob.glob(f"{args.summaries_dir}/*.summary")

    for summary_table in summary_tables[:10]:
        print(summary_table)
        add_lineage_to_table(summary_table)







if __name__ == "__main__":
    main()
