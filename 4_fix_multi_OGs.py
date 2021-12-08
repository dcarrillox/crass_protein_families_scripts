'''
This script reads the .tsv file where I classified as 'merge' or 'leave' the
different proteins ("3_multi_OGs_single-copy_hmmseach_reviewed.tsv"), and returns
the OG_xx.faa file after "leaving" or "merging" the proteins.

Lastly, it returns a flat text file 'all_leave_OGs.txt' with the OGs where only
leave proteins were detected.
'''

import argparse
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--classified_tsv_file',
                               dest='tsv_file',
                               required=True,
                               help='"3_multi_OGs_single-copy_hmmseach_reviewed.tsv" file '
                               'with the "leave" or "merge" classifications'
                               )

    requiredArgs.add_argument('-gt', '--genome_tables_dir',
                               dest='genome_tables_dir',
                               required=True,
                               help='Directory with the genome tables.'
                               )

    requiredArgs.add_argument('-og', '--og_faa_dir',
                               dest='og_faa_dir',
                               required=True,
                               help='filtered OG .faa files'
                               )


    requiredArgs.add_argument('-o', '--out_dir',
                               dest='out_dir',
                               required=True,
                               help='directory to put the new .faa files'
                               )

    return parser.parse_args()


def tsv_to_dicts(tsv_file):
    '''
    Reads the tsv classification file and returns two dictionaries, one with the
    'merge' proteins and other with the 'leave' ones. k=OG_id, v=[proteins_id]
    '''
    # read tsv file and retain only 'merge' or 'leave' proteins
    lines = [line.replace("\n", "").split("\t") for line in open(tsv_file).readlines() if line.split("\t")[-3] in ["merge", "leave"]]
    
    # separate in to 'merge' and 'leave' dictionaries
    # initiate dicts
    merge = {line[0]:list() for line in lines if line[-3] == "merge"}
    leave = {line[0]:list() for line in lines if line[-3] == "leave"}

    # fill in the dicts. k=OG, v=[proteins_ids]
    for line in lines:
        if line[-3] == "merge":
            merge[line[0]].append(line[1:-3])
        else:
            leave[line[0]].append(line[1:-3])

    return merge, leave

def get_start_stop_strand(genome_tables_dir):
    '''
    Reads genome tables and stores the start, stop and strand info per protein
    '''
    protein_info = dict()

    # get genome tables
    tables = glob.glob(f"{genome_tables_dir}/*table")

    for table in tables:
        lines = [line.split("\t") for line in open(table).readlines()[1:]]
        # iterate proteins and store start, stop and strand
        for line in lines:
            protein_info[line[-2]] = [line[2], line[3], line[4]]


    # fix HvCF_A4_ms_1|93|44. Is negative strand, but I still want to merge with the other
    # part of the protein since diamond confirms it is the final part of the protein
    #protein_info["HvCF_A4_ms_1|93|44"][2] = "+"

    return protein_info

def merge_proteins_ids(merge_dict, protein_info):
    '''
    Iterate the proteins 'merge'
    '''

    to_merge = dict()

    for og, og_proteins in merge_dict.items():
        to_merge[og] = list()
        # separate proteins into the genomes they come from
        # first initiate a dict with the genomes id
        genomes = {genome:list() for genome in list(set([protein[0].split("|")[0] for protein in og_proteins]))}

        # now separate by genome
        for protein in og_proteins:
            # store only the protein identifier
            genomes[protein[0].split("|")[0]].append(protein[0])

        # now process the proteins in each genome
        for genome, genome_proteins in genomes.items():
            # check that the strands of all the proteins are the same
            strands = [protein_info[protein][2] for protein in genome_proteins]
            strand  = list(set(strands))

            if len(strand) == 1:
                # remove redundant proteins product of multidomain hit in the protein
                uniq = list(set(genome_proteins))
                # now sort the proteins according to the strand:
                # negative, late genes go first; positive, earlier genes go first
                if strand[0] == "-":
                    sorted_proteins = sorted(uniq, key=lambda protein: int(protein.split("|")[-1]), reverse=True)
                else:
                    sorted_proteins = sorted(uniq, key=lambda protein: int(protein.split("|")[-1]))

                # the script assumes that there is only one split protein to merge,
                # ie. for that OG there is only one "group" of ORFs to merge in the genome.
                # if the distance between ORFs is larger than 3 genes I assume that
                # it is an ORF split at the ends of the genome and the order of merging
                # has to be inverted.
                invert = False
                for i in range(0, len(sorted_proteins)- 1):
                    current = int(sorted_proteins[i].split("|")[-1])
                    next    = int(sorted_proteins[i+1].split("|")[-1])

                    if abs(current - next) > 40:
                        invert = True
                        print(sorted_proteins)

                # if the distance between genes is not that large, merge them according to the order in the list,
                # that was given by the strand

                # I do the "OBAG01000118|121|108" gene manually since it is split and
                # at the ends of the genome.
                if "OBAG01000118|121|108" in sorted_proteins:
                    to_merge[og].append(["OBAG01000118|165|4", "OBAG01000118|153|2", "OBAG01000118|121|108"])

                else:
                    if not invert:
                        to_merge[og].append(sorted_proteins)
                    else:
                        sorted_proteins.reverse()
                        to_merge[og].append(sorted_proteins)

            else:
                print(strands, genome_proteins)

    return to_merge



def main():

    args = parse_args()

    # separate into 'merge' and 'leave' dictionaries
    merge, leave = tsv_to_dicts(args.tsv_file)

    # read genome tables to get start, stop and strand information of the proteins
    # dictionary, k=protein_id, v=[start, stop, strand]
    protein_info = get_start_stop_strand(args.genome_tables_dir)

    # figure out how the proteins should be merged
    to_merge = merge_proteins_ids(merge, protein_info)
    # write a report file with merging results ids
    with open(f"{args.out_dir}/merging_summary.txt", "w") as fout:
        for og, genomes in to_merge.items():
            for proteins in genomes:
                fout.write(f"{og}\t" + ", ".join(proteins) + "\n")




    # iterate the OGs while opening its filtered faa file and merging sequences
    for og, genomes in to_merge.items():
        # parse OG filtered faa file, store to dict: k=protein_id, v=record
        records = SeqIO.parse(f"{args.og_faa_dir}/{og.replace('_singlec','')}.faa", "fasta")
        proteins_seqs = {record.id:record for record in records}

        remove = list()
        for proteins in genomes:
            genome_id = proteins[0].split("|")[0]

            # merge the n_gene of the proteins
            joint_n = "_".join([protein.split("|")[-1] for protein in proteins])
            seq = str()
            for protein in proteins:
                seq += str(proteins_seqs[protein].seq).replace("*", "")
                remove.append(protein)

            # create the protein_id for the joint protein_
            joint_id = f"{genome_id}|{len(seq)}|{joint_n}"
            new_record = SeqRecord(Seq(seq), id=joint_id, description="")

            proteins_seqs[joint_id] = new_record

        tow = list()
        for protein, record in proteins_seqs.items():
            if protein not in remove:
                tow.append(record)

        with open(f"{args.out_dir}/{og.replace('_singlec','')}.faa", "w") as fout:
            SeqIO.write(tow, fout, "fasta")



    # write the flat file 'all_leave_OGs.txt'.
    # iterate the OGs in 'leave' and check if they are also in 'merge'. If not, write it to the file.
    with open(f"{args.out_dir}/all_leave_OGs.txt", "w") as fout:
        for OG, proteins in leave.items():
            if OG not in merge:
                fout.write(f"{OG}\n")


if __name__ == "__main__":
    main()
