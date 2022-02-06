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

    # create dir
    encode_dir = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{args.run_name}/3_vs2/1_vs2_results"
    os.makedirs(encode_dir, exist_ok=True)

    # copy results to encode
    os.system(f"cp final-viral-boundary.tsv {encode_dir}/final-viral-boundary_0{args.part_name}.tsv")
    os.system(f"cp final-viral-score.tsv {encode_dir}/final-viral-score_0{args.part_name}.tsv")


if __name__ == "__main__":
    main()
