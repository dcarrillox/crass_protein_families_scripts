'''
Process the summary tables (one per CL) to create a table for the rank:

- superkingdom
- phylum
- class
- order
- family
- genus

In these tables, rows are the CLs, columns are the different taxas. For each of
taxas, the value represents its fraction among all the TRUE hits in the CL.

It first concat all the summaries in one df, from which it gets all the rank names.
Then, for each of the ranks, it creates a df with those names and processes the tables
'''

import pandas as pd
import numpy as np
import os, glob, argparse, multiprocessing, time, itertools
from functools import partial



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-s', '--summary_dir',
                               dest='summary_dir',
                               required=True,
                               help='directory to save the .summary files'
                               )
    requiredArgs.add_argument('-o', '--outdir',
                               dest='outdir',
                               required=True,
                               help='directory to save the .faa files'
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

def read_summary_table(table_file):
    '''
    Reads the table file and keeps only included (TRUE) hits
    '''
    df = pd.read_csv(table_file, sep="\t", header=0, index_col=0, low_memory=False)
    df = df[df.included]
    return df

def count_names_in_summary_table(table_file, rank, rank_df_counts, rank_df_perctg, crass_ranks_names):
    cl_id = os.path.basename(table_file).split(".")[0]
    df = pd.read_csv(table_file, sep="\t", header=0, index_col=0, low_memory=False)
    df = df[df.included]
    if not df.empty:
        nhits = df.included.sum()

        df[rank] = df[rank].fillna("Unclassified")

        # count number of ncbi_crass
        ncbi_crass_counts = df.ncbi_crass.sum()
        rank_df_counts.loc[cl_id, "ncbi_crass"] = ncbi_crass_counts
        rank_df_perctg.loc[cl_id, "ncbi_crass"] = round(ncbi_crass_counts/nhits, 3) * 100

        # iterate the the counts for each name and fill in the rank_df
        for name, count in df[rank].value_counts().iteritems():
            rank_df_counts.loc[cl_id, name] = count
            rank_df_perctg.loc[cl_id, name] = round(count/nhits, 3) * 100

        # subtract ncbi_crass if Viruses, Uroviricota, Caudoviricetes, Caudovirales, Podoviridae
        if rank != "genus":
            crass_name = crass_ranks_names[rank]
            rank_df_counts.loc[cl_id, crass_name] = rank_df_counts.loc[cl_id, crass_name] - ncbi_crass_counts
            rank_df_perctg.loc[cl_id, name] = round(rank_df_counts.loc[cl_id, crass_name]/nhits, 3) * 100



def main():

    args = parse_args()


    summary_files = glob.glob(f"{args.summary_dir}/*.summary")
    cl_ids = sorted([os.path.basename(file).split(".")[0] for file in summary_files])

    print("parsing summary files...")
    start = time.time()
    # process hmmsearch files in parallel
    # cls_seqsids is a list of lists, in each list the first item is the cl_id, the rest are hits protein ids
    pool = multiprocessing.Pool(processes=18)
    all_summaries_dfs = pool.map(read_summary_table, summary_files)
    pool.close()
    pool.join()
    print("Done")
    timeit(start)
    print()

    print("concat in one df to get all the possible taxa names...")
    start = time.time()
    all_summaries = pd.concat(all_summaries_dfs, axis=0, ignore_index=False)
    all_summaries = all_summaries.fillna("Unclassified")
    timeit(start)
    print()

    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus"]

    # create a dict with the taxa names for crass, from superkingdom to family.
    crass_names = ["Viruses", "Uroviricota", "Caudoviricetes", "Caudovirales", "Podoviridae"]
    crass_ranks_names = {rank:name for rank, name in zip(target_ranks[:-1], crass_names)}


    for rank in target_ranks:
        print(f"Generating {rank} table...")
        start = time.time()
        # get the unique names for the rank
        names = sorted(list(all_summaries[rank].unique()) + ["ncbi_crass"])
        rank_df_counts = pd.DataFrame(0, columns=names, index=cl_ids)
        rank_df_perctg = pd.DataFrame(0, columns=names, index=cl_ids)

        for file in summary_files:
            count_names_in_summary_table(file, rank, rank_df_counts, rank_df_perctg, crass_ranks_names)

        # rank_df_counts = rank_df_counts.T
        # rank_df_perctg = rank_df_perctg.T

        rank_df_counts.to_csv(f"{args.outdir}/{rank}_count.txt", index=True, header=True, sep="\t")
        rank_df_perctg.to_csv(f"{args.outdir}/{rank}_perctg.txt", index=True, header=True, sep="\t")

        timeit(start)
        print()

if __name__ == "__main__":
    main()
