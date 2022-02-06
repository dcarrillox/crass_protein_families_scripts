'''

'''

import os, glob, argparse, multiprocessing, time
from functools import partial
import pandas as pd
import numpy as np
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



def get_CLs_to_analyze():
    '''
    Get cl_ids from MRCAs directories
    '''

    mrcas_files = glob.glob("/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/0_mrcas/*.txt")
    cl_ids = [os.path.basename(file).replace(".txt", "") for file in mrcas_files]

    return cl_ids

def get_contigs_status():
    '''
    '''
    #
    # # list .contigs files, get identifiers from them and init a dict
    contigs_files = glob.glob("/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/1_ipg/*.contigs")
    contigs_ids = list()
    for file in contigs_files:
        contigs_ids += [line.strip() for line in open(file).readlines()]
    # init them as "not_downloaded"
    contigs_status = {contig:"not_down" for contig in contigs_ids}


    # get nr length stats from different rounds
    contigs_stats_files = glob.glob("/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/3_vs2/0_fasta/nr_contigs.stats")

    #contigs_status = dict()
    for stats_file in contigs_stats_files:
        lines = [line.strip().split("\t") for line in open(stats_file).readlines()[1:]]
        # split by space the first column with the header to have only the contig_id
        for line in lines:
            contig_id = line[0].split(" ")[0]
            length = int(line[1])

            if length < 1500:
                contigs_status[contig_id] = "short"
            else:
                contigs_status[contig_id] = "ok"


    return contigs_status

def read_summary_file(summary_table):
    '''
    Read summary file to know for which proteins should I create an entry in the final table
    '''

    summary_df = pd.read_csv(summary_table, sep="\t", header=0, index_col=0, low_memory=False)
    # keep only Bacteria
    bacteria_summary_df = summary_df[(summary_df["included"] == True) & (summary_df["superkingdom"] != "Viruses")]
    prots_ids = bacteria_summary_df.index.to_list()

    return prots_ids

def get_ipgs_for_prots(cl_id, prots_ids):
    '''
    Checks the .ipg directories of all analysis rounds, concat all of them in the same df (or list of lists).
    Then, parse the merged ipg table and assign the unique identifiers associated with each protein
    '''

    # first get which is the ipg for the protein that is in the tree. NB this does
    # not need to be the one in the .ipg_filtered. Because of that I need to first


    # list all the ipg tables of different rounds for the cl
    ipg_tables = glob.glob(f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/1_ipg/{cl_id}.ipg") # .ipg_filtered

    ipg_lines = list()
    for ipg_table in ipg_tables:
        lines = [line.strip().split("\t") for line in open(ipg_table, "r").readlines()[1:]]
        ipg_lines += lines

    # iterate the ipg table(s) looking for the unique identifiers for each protein
    prots_ipgs = dict()
    for line in ipg_lines:
        if line[6] in prots_ids:
            prots_ipgs[line[6]] = line[0]
            if line[6] == "MCC6012565.1":
                print("whut")

    # do the inverse, store in a dict k=ipg_id v=prot_id
    ipgs_prots = dict()
    for prot, ipg_id in prots_ipgs.items():
        ipgs_prots[ipg_id] = prot

    return prots_ipgs, ipgs_prots

def get_contigs_ipgs(cl_id, ipgs_prots):
    '''
    For each ipg, store the contigs associated with it.
    '''

    # list all the ipg tables of different rounds for the cl
    ipg_tables = glob.glob(f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/1_ipg/{cl_id}.ipg_filtered")

    ipg_lines = list()
    for ipg_table in ipg_tables:
        lines = [line.strip().split("\t") for line in open(ipg_table, "r").readlines()[1:]]
        ipg_lines += lines

    # iterate the .ipg lines while storing the contigs downloaded for each bacteria protein
    ipgs_contigs = {ipg_id:dict() for ipg_id in ipgs_prots}
    for line in ipg_lines:
        # cl_0465-MCC6012565.1 check
        if line[0] in ipgs_contigs:
            if line[2] in ipgs_contigs[line[0]]:
                ipgs_contigs[line[0]][line[2]] += [[int(line[3]), int(line[4])]]
            else:
                ipgs_contigs[line[0]][line[2]] = [int(line[3]), int(line[4])]
        else:
            print(f"{line[0]}, {cl_id} is not in the contigs dictionary")

    return ipgs_contigs

def parse_vs2_results():
    '''
    List all the vs2 results from different rounds and store them in a dict
    '''

    vs2_dirs = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/3_vs2/1_vs2_results"

    # read scores files
    vs2_scores_files = glob.glob(f"{vs2_dirs}/final-viral-score*.tsv")

    vs2_scores_lines = list()
    for file in vs2_scores_files:
        lines = [line.strip().split("\t") for line in open(file, "r").readlines()[1:]] # discard header
        vs2_scores_lines += lines

    # store contigs in the results in a dict along with its score
    vs2_results_dict = {line[0].split("|")[0]:dict() for line in vs2_scores_lines}
    for line in vs2_scores_lines:
        contig = line[0].split("|")[0]
        vs2_results_dict[contig][line[0]] = [float(line[3])]
        # control for lt2gene contigs: set them as viral score 1, coordinates 1 to 10000000
        if "||lt2gene" in line[0]:
            vs2_results_dict[contig][line[0]] = [1, 1, 10000000]


    # read boundaries files
    vs2_boundary_files = glob.glob(f"{vs2_dirs}/final-viral-boundary*.tsv")

    vs2_boundary_lines = list()
    for file in vs2_boundary_files:
        lines = [line.strip().split("\t") for line in open(file, "r").readlines()[1:]] # discard header
        vs2_boundary_lines += lines

    # store start and stop prophage coordinates in each contig
    for line in vs2_boundary_lines:
        contig = line[0]
        # not all the contigs in boundary.tsv are in score.tsv. Weird, don't know why yet
        # check if that is the case
        if contig in vs2_results_dict:
            vs2_results_dict[contig][line[-1]] += [int(line[3]), int(line[4])]
        else:
            print(f"Contig {contig} not in scores file.")
            vs2_results_dict[contig] = {line[-1]:[0, int(line[3]), int(line[4])]}


    return vs2_results_dict

def check_protein_in_prophage(prot_start, prot_end, prophage_start, prophage_end):
    if prophage_start < prot_start < prophage_end or prophage_start < prot_end < prophage_end:
        return True
    else:
        return False

def classify_ipg(cl_id, bact_proteins, ipgs_prots, ipgs_contigs, vs2_results_dict, contigs_status):
    '''
    '''

    # init a dataframe with the bacteria proteins as index
    columns = ["ipg_id", "prophage", "discard", "n_contigs", "n_prophage", "n_short", "n_notdown", "contigs", "prophage_contigs", "short_contigs", "notdown_contigs"]
    final_df = pd.DataFrame(0, index=bact_proteins, columns=columns)
    final_df["prophage"] = False
    final_df["discard"] = False
    final_df["contigs"] = ""
    final_df["short_contigs"] = ""
    final_df["prophage_contigs"] = ""


    # iterate the ipgs
    for ipg_id, prot_id in ipgs_prots.items():
        final_df.loc[prot_id, "ipg_id"] = ipg_id

        # list to store the prophage contigs ids, so I can write them later together
        prophage_contigs = list()

        # check if there are downloaded contigs associated with the ipg_id
        if ipgs_contigs[ipg_id]:
            # set the number of contigs inspected
            final_df.loc[prot_id, "n_contigs"] = len(ipgs_contigs[ipg_id])
            final_df.loc[prot_id, "contigs"] = ",".join(ipgs_contigs[ipg_id])

            # get short contigs number and ids
            short_contigs = [contig_id for contig_id in ipgs_contigs[ipg_id] if contigs_status[contig_id] == "short"]
            final_df.loc[prot_id, "n_short"] = len(short_contigs)
            final_df.loc[prot_id, "short_contigs"] = ",".join(short_contigs)

            # get not_downloaded contigs number and ids
            notdown_contigs = [contig_id for contig_id in ipgs_contigs[ipg_id] if contigs_status[contig_id] == "not_down"]
            final_df.loc[prot_id, "n_notdown"] = len(notdown_contigs)
            final_df.loc[prot_id, "notdown_contigs"] = ",".join(notdown_contigs)

            # iterate the contigs associated with each ipg
            for contig_id, prot_coords in ipgs_contigs[ipg_id].items():

                # first check it is not a "short" contig
                if contigs_status[contig_id] != "short":
                    check = False
                    if contig_id in vs2_results_dict:
                        #print(cl_id, ipg_id, contig_id)
                        for prophage in vs2_results_dict[contig_id]:
                            #print(prophage, vs2_results_dict[contig_id][prophage])
                            score = vs2_results_dict[contig_id][prophage][0]
                            prophage_start = vs2_results_dict[contig_id][prophage][1]
                            prophage_end = vs2_results_dict[contig_id][prophage][2]

                            if check_protein_in_prophage(prot_coords[0], prot_coords[1], prophage_start, prophage_end):
                                check = True
                        if check:
                            final_df.loc[prot_id, "n_prophage"] += 1
                            final_df.loc[prot_id, "prophage"] = True
                            prophage_contigs.append(contig_id)

            # write prophage contigs, if any
            final_df.loc[prot_id, "prophage_contigs"] = ",".join(prophage_contigs)


        # Protein was not inspected in this round
        else:
            final_df.loc[prot_id, :] = np.nan

    return final_df

def process_cl(contigs_status, vs2_results_dict, run_name, cl_id):
    '''
    input: .ipg table

    1. Read summary table
    2. Read all .ipg tables for the CL from different rounds
    '''

    print(cl_id)

    # parse summary table, get the list of proteins in the tree
    summary_dir = "/home/danielc/projects/Bas_phages/5_nr_screening/2_hits_summary"
    summary_table = f"{summary_dir}/{cl_id}.summary"
    bact_proteins = read_summary_file(summary_table)
    if "MCC6012565.1" in bact_proteins:
        print("OK")

    # for each of the proteins, get its ipg_id, and for each ipg_id get its proteins
    prots_ipgs, ipgs_prots = get_ipgs_for_prots(cl_id, bact_proteins)

    # associate each ipg with the contigs where they have been found
    ipgs_contigs = get_contigs_ipgs(cl_id, ipgs_prots)

    # iterate the ipg ids while assessing if they come from a prophage
    final_df = classify_ipg(cl_id, bact_proteins, ipgs_prots, ipgs_contigs, vs2_results_dict, contigs_status)


    # write df to file
    outfile = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{run_name}/3_vs2/2_classified/{cl_id}.prophages"
    final_df.to_csv(outfile, sep="\t", header=True, index=True)





def main():

    args = parse_args()

    # get the CL_ids to analyze from the 2_trees folder
    cl_ids = get_CLs_to_analyze()

    # parse contigs .fasta files to know which ones are short or could not be downloaded
    print("Parsing contigs .fasta files...")
    contigs_status = get_contigs_status()

    # store all the vs2 results I have so far from different rounds in a dict
    print("Parsing vs2 results...")
    vs2_results_dict = parse_vs2_results()


    # part = partial(process_cl, contigs_status, vs2_results_dict, args.run_name)
    # pool = multiprocessing.Pool(processes=18)
    # pool.map(part, cl_ids)
    # pool.close()
    # pool.join()

    for cl_id in cl_ids:
        process_cl(contigs_status, vs2_results_dict, args.run_name, cl_id)

if __name__ == "__main__":
    main()
