
import os, glob, argparse, multiprocessing, time
from functools import partial
from Bio import Entrez
import pandas as pd
import numpy as np
from math import ceil
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-n', '--run_name',
                               dest='run_name',
                               required=True,
                               help='name of the analysis (round1, round2...)'
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    # create oudirs
    root_outdir = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial"
    os.makedirs(f"{root_outdir}/{args.run_name}/3_vs2/0_fasta", exist_ok=True)

    # concat all contigs .fasta files in one file
    print(f"Concatenating .fasta files to '{root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta'...")
    fasta_files = glob.glob(f"{root_outdir}/{args.run_name}/2_fasta/*.fasta")

    to_write = list()
    for fasta_file in fasta_files:
        records = [record for record in SeqIO.parse(fasta_file, "fasta")]
        to_write += records

    with open(f"{root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta", "w") as fout:
        SeqIO.write(to_write, fout, "fasta")

    #print(f"Concatenating .fasta files to '{root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta'...")
    #print(f"cat {' '.join(fasta_files)} > {root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta")
    #os.system(f"cat {' '.join(fasta_files)} > {root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta")

    # remove redundancy
    print("Removing redundancy...")
    os.system(f"seqkit rmdup -j 18 < {root_outdir}/{args.run_name}/3_vs2/0_fasta/all_contigs.fasta > {root_outdir}/{args.run_name}/3_vs2/0_fasta/nr_contigs.fasta")

    # get length stats
    print("Calculating lengths...")
    os.system(f"seqkit fx2tab --length --name --header-line < {root_outdir}/{args.run_name}/3_vs2/0_fasta/nr_contigs.fasta > {root_outdir}/{args.run_name}/3_vs2/0_fasta/nr_contigs.stats")

    # split in 6 parts
    print("Splitting file...")
    os.system(f"seqkit split2 -j10 -p 10 {root_outdir}/{args.run_name}/3_vs2/0_fasta/nr_contigs.fasta")


if __name__ == "__main__":
    main()
