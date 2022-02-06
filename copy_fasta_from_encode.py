import os, glob, multiprocessing, time, itertools
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-n', '--run_name',
                               dest='run_name',
                               required=True,
                               help='name of the analysis (round1, round2...)'
                               )
    requiredArgs.add_argument('-p', '--part_name',
                               dest='part_name',
                               required=True,
                               help='just the number: 1, 2, 3...'
                               )


    return parser.parse_args()


def main():

    args = parse_args()

    # remove previous vs2 run
    os.system("rm -rf config.yaml *.tsv *.fasta iter-0 log final-viral-combined.fa")

    # copy .fasta part file
    part_file = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{args.run_name}/3_vs2/0_fasta/nr_contigs.fasta.split/nr_contigs.part_00{args.part_name}.fasta"
    os.system(f"cp {part_file} .")










if __name__ == "__main__":
    main()
