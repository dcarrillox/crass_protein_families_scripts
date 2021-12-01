'''
Parses the .tab files generated with parse_hhblits_tab.py to create
the input, three-columns tabular file, needed for MCL.
'''

import glob, os, argparse
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hhr_tab_dir',
                               dest='hhr_tab_dir',
                               required=True,
                               help='directory with the .hhsearch files'
                               )
    requiredArgs.add_argument('-m', '--msa_dir',
                               dest='msa_dir',
                               required=False,
                               help='directory with the .mafft files, used '
                               'to make get the relation sequence-PF_id.'
                               )
    requiredArgs.add_argument('-omcl', '--omcl',
                               dest='omcl',
                               required=False,
                               help='output 3-tab file for MCL.'
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # list hht_tab files
    hhr_tab_files = glob.glob(f"{args.hhr_tab_dir}/*.tab")

    print("Reading .tab files...")
    df_list = [pd.read_csv(hhr_tab_file, sep="\t", header=0) for hhr_tab_file in hhr_tab_files]
    #concatenate them together
    print("Done. Concatenating ")
    hhr_tab_df = pd.concat(df_list, ignore_index=True)
    # remove probs <95% and qcov <0.5
    hhr_tab_df = hhr_tab_df[(hhr_tab_df.prob.ge(95)) &
                            (hhr_tab_df.qcov.ge(0.5)) &
                            (hhr_tab_df.tcov.ge(0.5)) &
                            (hhr_tab_df["query"] != hhr_tab_df["target"])]

    #
    to_write = hhr_tab_df[["query", "target", "qscore"]]
    to_write.to_csv(args.omcl, sep="\t", header=None, index=False)



if __name__ == "__main__":
    main()
