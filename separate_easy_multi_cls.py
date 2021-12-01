'''
Given the .faa files for the cls obtained with create_cls_faa.py, separates the cls
in two folders: multi and single. For the multi ones, it also creates the single version
of them and a .faa file with the multi proteins, to align the single version and compare
to the multi proteins.
'''

import argparse, glob, os
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--in_faas_dir',
                               dest='in_faas_dir',
                               required=True,
                               help='dir with the cl_XX.faa files ("Bas_phages/4_protein_families/0_faa")'
                               )
    requiredArgs.add_argument('-osd', '--out_single_dir',
                               dest='out_single_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-omd', '--out_multi_dir',
                               dest='out_multi_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-t', '--out_tab_file',
                               dest='out_tab_file',
                               required=True,
                               help='out 2-columns file saying which cls are single or multi'
                               )

    return parser.parse_args()


def main():

    args = parse_args()

    # list the .faa files
    cls_faas = glob.glob(f"{args.in_faas_dir}/*.faa")

    # create the output directories
    os.makedirs(args.out_multi_dir, exist_ok=True)
    os.makedirs(f"{args.out_multi_dir}/0_singlec", exist_ok=True)
    os.makedirs(args.out_single_dir, exist_ok=True)
    os.makedirs(f"{args.out_single_dir}/0_faa", exist_ok=True)

    # create the list to store the separate_summary
    summary = [["CL", "type"]]

    # create the list to store all the multi proteins
    multi_to_write = list()

    # iterate the cls .faas and process them
    for cl_faa in cls_faas:
        cl_id = os.path.basename(cl_faa).split(".faa")[0]
        # parse to dict
        records = SeqIO.to_dict(SeqIO.parse(cl_faa, "fasta"))

        # get the genomes in the cl
        genomes_cl = [record.split("|")[0] for record in records]
        uniq_genomes_cl = set(genomes_cl)

        # check if duplicates
        # non duplicated genomes = single
        if len(genomes_cl) == len(uniq_genomes_cl):
            os.system(f"cp {cl_faa} {args.out_single_dir}/0_faa")
            summary.append([cl_id, "single"])
        else:
            summary.append([cl_id, "multi"])
            # identify the duplicated genomes
            dupl = [genome for genome in uniq_genomes_cl if genomes_cl.count(genome) > 1]

            single_version = [record for record_id, record in records.items() if record_id.split("|")[0] not in dupl]
            multi_proteins = [record for record_id, record in records.items() if record_id.split("|")[0] in dupl]
            multi_to_write += multi_proteins

            with open(f"{args.out_multi_dir}/0_singlec/{cl_id}_singlec.faa", "w") as fout:
                SeqIO.write(single_version, fout, "fasta")

    # write to .faa all the multi_proteins
    with open(f"{args.out_multi_dir}/proteins_multi.faa", "w") as fout:
        SeqIO.write(multi_to_write, fout, "fasta")

    # write the summary file
    with open(args.out_tab_file, "w") as fout:
        for line in summary:
            fout.write("\t".join(line) + "\n")




if __name__ == "__main__":
    main()
