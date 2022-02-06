'''
It merges the filtered .faa files obtained from NR with extract_seqs_NR.py & filter_nr_hits.py
with the .faa sequences of my CLs.
'''

import pandas as pd
import numpy as np
import os, glob, argparse, multiprocessing, time, itertools
from functools import partial



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-nf', '--ncbi_faa_dir',
                               dest='ncbi_faa_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-cf', '--cls_faa_dir',
                               dest='cls_faa_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-s', '--summaries_dir',
                               dest='summaries_dir',
                               required=True,
                               help=''
                               )

    requiredArgs.add_argument('-o', '--outdir',
                               dest='outdir',
                               required=True,
                               help='directory to save the .faa files'
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # check which CLs got Bacteria hits
    bacteria_cls = list()
    for table in glob.glob(f"{args.summaries_dir}/*.summary"):
        df = pd.read_csv(table, header=0, sep="\t", index_col=0, low_memory=False)
        df.fillna("", inplace=True)
        # keep only included
        included   = df.loc[df.included, :]
        n_bacteria = included.loc[df.superkingdom.isin(["Bacteria", ""]), 'superkingdom'].count()

        if n_bacteria > 0:
            cl_id = os.path.basename(table).split(".")[0]
            bacteria_cls.append(cl_id)


    print(f"{len(bacteria_cls)} CLs got Bacteria hits")
    # iterate the CLs with Bacteria hits
    for cl in bacteria_cls:
        cl_faa_file   = f"{args.cls_faa_dir}/{cl}.faa"
        ncbi_faa_file = f"{args.ncbi_faa_dir}/{cl}.faa"
        outfile       = f"{args.outdir}/{cl}_ncbi.faa"

        os.system(f"cat {cl_faa_file} {ncbi_faa_file} > {outfile}")


    # # check which CLs have been already merged
    # done = [os.path.basename(cl).split("_ncbi")[0] for cl in glob.glob(f"{args.outdir}/*_ncbi.faa")]
    #
    # # list faa files obtained from NR
    # ncbi_faa_files = glob.glob(f"{args.ncbi_faa_dir}/*.faa")
    # print(f"{len(ncbi_faa_files)} in the ncbi folder")
    # # remove those already processed
    # ncbi_faa_files = [faa_file for faa_file in ncbi_faa_files if os.path.basename(faa_file).split(".")[0] not in done]
    # print(f"{len(ncbi_faa_files)} will be processed")
    #
    # print("Merging NCBI and Crassvirales .faa files...")
    # for ncbi_faa_file in ncbi_faa_files:
    #     cl_id = os.path.basename(ncbi_faa_file).split(".")[0]
    #     crassvirales_faa_file = f"{args.cls_faa_dir}/{cl_id}.faa"
    #     outfile = f"{args.outdir}/{cl_id}_ncbi.faa"
    #     os.system(f"cat {crassvirales_faa_file} {ncbi_faa_file} > {outfile}")
    #
    #
    #

if __name__ == "__main__":
    main()
