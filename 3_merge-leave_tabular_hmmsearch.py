'''
I don't want to process by hand the table from 2_hmmer3-text_to_tabular, again.
It returns the same table but with the merge/leave labels added.
'''


import argparse
import glob
import multiprocessing
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_tabular_file',
                               dest='hmmsearch_tabular_file',
                               required=True,
                               help=' ".parsed-oghit" file from "2_hmmer3-text_to_tabular.py" '
                               )

    requiredArgs.add_argument('-gt', '--genome_tables_dir',
                               dest='genome_tables_dir',
                               required=True,
                               help=''
                               )


    requiredArgs.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )

    return parser.parse_args()


def parse_genome_tables_strand(genome_tables_dir):
    '''
    '''
    genome_tables = glob.glob(f"{genome_tables_dir}/*.table")

    genes_strand = dict()

    for table in genome_tables:
        lines = [line.split("\t") for line in open(table).readlines()[1:]]

        for line in lines:
            genes_strand[line[-2]] = line[4]

    return genes_strand

def check_protein_presence(genome_hits):
    '''
    Checks if a protein is present more than 1 time in the hits. If so, don't process
    the genome.
    '''
    # get the unique protein ids
    uniq_prots_ids = list(set([hit[1] for hit in genome_hits]))
    # check the length
    if len(uniq_prots_ids) == len(genome_hits):
        return True
    else:
        return False

def process_og_hits(ogs_hits, genes_strand, og_id):
    '''
    Function to process the lines for a given OG.
    '''

    # group the hits per genome. First get the unique genome_ids
    #print(ogs_hits[og_id])
    genomes_ids = list(set([line[1].split("|")[0] for line in ogs_hits[og_id]]))
    og_genomes = {genome_id:list() for genome_id in genomes_ids}

    for hit in ogs_hits[og_id]:
        genome_id = hit[1].split("|")[0]
        # remove the "\n" from the line. I didn't this before to keep the length of the lines
        hit[-1] = hit[-1].replace("\n", "")
        og_genomes[genome_id].append(hit)


    # process only the 2 hits genomes
    for genome, hits in og_genomes.items():
        #print(hits)
        if check_protein_presence(hits) and len(hits) == 2:
            if hits[0][2] not in ["No hit", "Other hits"] and hits[1][2] not in ["No hit", "Other hits"]:
                # check that genes are separated by less than 3 genes
                n1 = hits[0][1].split("|")[-1]
                n2 = hits[1][1].split("|")[-1]
                if int(n2) - int(n1) < 3:

                    # check strands
                    strands = list(set([genes_strand[hits[0][1]], genes_strand[hits[1][1]]]))
                    if len(strands) == 1:
                        # check if the align over the whole profile
                        profile_len = int(hits[0][4])
                        cov_1 = float((int(hits[0][8]) - int(hits[0][7])) / profile_len)
                        cov_2 = float((int(hits[1][8]) - int(hits[1][7])) / profile_len)
                        # if they do, add "leave" and "whole_profile" to the hits
                        if cov_1 >= 0.7 and cov_2 >= 0.7:
                            hits[0] += ["leave", "whole_profile"]
                            hits[1] += ["leave", "whole_profile"]

                        # if they don't...
                        else:
                            # if strand is "-" revert the order of the proteins
                            if strands[0] == "-":
                                hits.reverse()

                            # check if they overlap
                            start_1 = int(hits[0][7])
                            end_1   = int(hits[0][8])
                            start_2 = int(hits[1][7])
                            end_2   = int(hits[1][8])

                            # if they don't overlap...
                            if end_1 < start_2 and cov_1+cov_2 >= 0.7:
                                hits[0] += ["merge", "non-overlap"]
                                hits[1] += ["merge", "non-overlap"]

                            # if they overlap...
                            elif end_1 > start_2 and cov_1+cov_2 >= 0.7:
                                # I merge them automatic if the overlap is < 10% of the profile length
                                overlap = float((end_1 - start_2) / profile_len)
                                if overlap <= 0.1:
                                    hits[0] += ["merge", "overlap"]
                                    hits[1] += ["merge", "overlap"]
                                else:
                                    continue
                                    #print(hits)

                    else:
                        #print(f"strand incongruence, {og_id} {genome}")
                        hits[0] += ["leave", "strand"]
                        hits[1] += ["leave", "strand"]


                elif 3 <= int(n2) - int(n1) < 5:
                    profile_len = int(hits[0][4])
                    cov_1 = float((int(hits[0][8]) - int(hits[0][7])) / profile_len)
                    cov_2 = float((int(hits[1][8]) - int(hits[1][7])) / profile_len)
                    if cov_1 >= 0.7 and cov_2 >= 0.7:
                        hits[0] += ["leave", "whole_profile"]
                        hits[1] += ["leave", "whole_profile"]

                elif int(n2) - int(n1) > 5:
                    hits[0] += ["leave", "distance"]
                    hits[1] += ["leave", "distance"]

            else:
                continue
                # to check



    to_return = list()
    for genome, hits in og_genomes.items():
        for hit in hits:
            to_return.append(hit)

    return to_return






def main():

    args = parse_args()

    # parse the genome tables to know in which strand is each gene
    genes_strand = parse_genome_tables_strand(args.genome_tables_dir)


    # read the tabular file
    lines = [line.split("\t") for line in open(args.hmmsearch_tabular_file).readlines()[1:]]

    # for each OG, assign the lines in a dict. k=genome  v=[prot1, prot2, prot3...]
    # first get the OG identifiers. Remove the ones for which I didn't make the alignment. Remove also the double empty lines in the tabular ("\n" after doing readlines())
    not_aligned = ["\n", "OG_131", "OG_161", "OG_2038", "OG_242", "OG_25", "OG_2602", "OG_264", "OG_2798", "OG_289", "OG_361", "OG_3814", "OG_4023", "OG_407", "OG_4247", "OG_4321", "OG_4349", "OG_4455", "OG_4599", "OG_4601", "OG_4624", "OG_471", "OG_4725", "OG_4726", "OG_4768", "OG_4799", "OG_4809", "OG_4810", "OG_4818", "OG_4827", "OG_4842", "OG_4868", "OG_488", "OG_499", "OG_503", "OG_526", "OG_63", "OG_65"]
    ogs_hits = {og:list() for og in list(set([line[0] for line in lines if line[0] not in not_aligned]))}
    # now assign the lines to the og. First remove the ogs above
    filt_lines = [line for line in lines if line[0] not in not_aligned]
    for line in filt_lines:
        ogs_hits[line[0]].append(line)


    partial_func = partial(process_og_hits, ogs_hits, genes_strand)
    pool = multiprocessing.Pool(processes=1)
    to_write = pool.map(partial_func, list(ogs_hits.keys()))
    #to_write = pool.map(partial_func, ["OG_1285_singlec"])
    pool.close()
    pool.join()



    # make the list flat
    to_write_flat = [line for og in to_write for line in og]
    # insert and "\n" between OGs to separate them in the out file
    last = to_write_flat[0][0]
    for line in to_write_flat:
        if line[0] != last:
            last = line[0]
            line[0] = "\n" + line[0]
        else:
            last = line[0]

    with open(args.out_file, "w") as fout:
        for line in to_write_flat:
            fout.write("\t".join(line) + "\n")







if __name__ == "__main__":
    main()
