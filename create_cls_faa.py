'''

'''

import argparse, glob, os
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--in_faa_file',
                               dest='in_faa_file',
                               required=True,
                               help='"2_Prodigal/3_final_annotation_formatted.faa"'
                               )
    requiredArgs.add_argument('-gt', '--genome_tables_dir',
                               dest='genome_tables_dir',
                               required=True,
                               help='to get the cl_(s_)id and protein_id from the genome table'
                               )
    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help='directory to put the CLs .faa files'
                               )
    return parser.parse_args()


def main():

    args = parse_args()

    # store all_proteins records in a dict
    all_prots_records = {record.id:record for record in SeqIO.parse(args.in_faa_file,"fasta")}

    # iterate genome tables, store proteins per CL in a dict, k=CL
    cls_prots = dict()
    for table in glob.glob(f"{args.genome_tables_dir}/*.table"):
        with open(table, "r") as fin:
            lines = [line.strip().split("\t") for line in fin.readlines()[1:]]
            for line in lines:
                # check the protein was assigned to a CL or SCL
                if "cl_" in line[-1]:
                    id = line[-1]
                    # store the protein along its CL or SCL
                    if id in cls_prots:
                        cls_prots[id] += [line[9]]
                    else:
                        cls_prots[id] = [line[9]]



    # iterate the CLs ids in the dictionary and write to .faa file
    for id, prots in cls_prots.items():
        to_write = [all_prots_records[prot] for prot in prots]
        outfile = f"{args.out_dir}/{id}.faa"
        with open(outfile, "w") as fout:
            SeqIO.write(to_write, fout, "fasta")



if __name__ == "__main__":
    main()
