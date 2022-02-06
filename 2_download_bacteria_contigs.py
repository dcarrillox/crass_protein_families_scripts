
import os, glob, argparse, multiprocessing, time
from functools import partial
from Bio import Entrez
import pandas as pd
import numpy as np
from math import ceil


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-n', '--run_name',
                               dest='run_name',
                               required=True,
                               help='name of the analysis (round1, round2...)'
                               )


    return parser.parse_args()



def get_already_downloaded_contigs():
    '''
    '''

    stats_files_dirs = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/3_vs2/0_fasta"
    stats_files = glob.glob(f"{stats_files_dirs}/nr_contigs.stats")

    stats_lines = list()
    for file in stats_files:
        lines = [line.strip().split("\t") for line in open(file).readlines()[1:]]
        stats_lines += lines


    # split first column by " " to get only the id of the contig
    already_down_contigs = set([line[0].split(" ")[0] for line in stats_lines])

    return already_down_contigs

def get_proteins_to_download(mrca_file):
    '''
    '''
    mrcas_df = pd.read_csv(mrca_file, sep="\t", header=0, index_col=0)
    mrcas_df = mrcas_df.fillna("NA")

    polytomies = [polytomy.split(",") for polytomy in mrcas_df["polytomies"].values.tolist()]
    sister_bacteria = [sister.split(",") for sister in mrcas_df["bacteria_seqs_sister"].values.tolist()]
    outgroups = [mrcas_df["outgroup"].tolist()]
    all = polytomies + sister_bacteria + outgroups
    flat_list = [item for sublist in all for item in sublist]
    bacteria_proteins = list(set(flat_list))

    if "NA" in bacteria_proteins:
        bacteria_proteins.remove("NA")

    return bacteria_proteins

def download_ipg(bacteria_proteins, cl_id, outdir, run_name):
    '''
    '''

    ipg_table = str()
    for i in range(0, len(bacteria_proteins), 200):
        batch = bacteria_proteins[i:i + 200]
        handle = Entrez.efetch(db="protein", id=batch, rettype="ipg", retmode="text")
        ipg_table += handle.read().decode("utf-8")
        handle.close()

    ipg_file = f"{outdir}/{run_name}/1_ipg/{cl_id}.ipg"
    with open(ipg_file, "w") as fout:
        if ipg_table:
            fout.write(ipg_table)

    return ipg_table

def filter_ipg(ipg_group):
    '''
    '''
    # change to list of lists
    ipg_group_list = ipg_group.values.tolist()

    # store which contigs I will consider for the download
    to_check = list()
    done = list()
    # first get RefSeq contigs
    for line in ipg_group_list:
        if line[1] == "RefSeq":
            to_check.append(line[2])
            done.append(line[2].replace("NZ_", ""))
	# process GenBank hits
    for line in ipg_group_list:
        if line[1] == "INSDC" and line[2] not in done:
            to_check.append(line[2])
            done.append(line[2])

    # filter df with the contigs in to_check
    to_check_boolean = ipg_group["Nucleotide Accession"].isin(to_check)
    ipg_group = ipg_group[to_check_boolean]

    # if there are more than 3 contigs, get the top 3 with the highest
    # start coordinate of the protein. With this I hope I am picking the longest
    # contigs.
    if ipg_group.shape[0] > 3:
        ipg_group["Start"] = pd.to_numeric(ipg_group["Start"])
        ipg_group = ipg_group.sort_values("Start", ascending=False)
        contigs_to_download = ipg_group.head(3)["Nucleotide Accession"].tolist()
        ipg_group = ipg_group.head(3)

    else:
        contigs_to_download = ipg_group["Nucleotide Accession"].tolist()

    return [ipg_group, contigs_to_download]

def parse_ipg(cl_id, outdir, run_name):
    '''
    '''

    ipg_file = f"{outdir}/{run_name}/1_ipg/{cl_id}.ipg"

    # read df, group by IPG identifier
    ipg_df = pd.read_csv(ipg_file, sep="\t", header=0)
    ipg_df = ipg_df[ipg_df["Nucleotide Accession"].notnull()]
    ipg_df = ipg_df.groupby("Id")

    parsed_ipg = [filter_ipg(group) for id, group in ipg_df]

    # get the filtered ipgs
    filtered_ipg_groups = [item[0] for item in parsed_ipg]
    filtered_ipg = pd.concat(filtered_ipg_groups)
    # write filtered ipg to file
    outfile = f"{outdir}/{run_name}/1_ipg/{cl_id}.ipg_filtered"
    filtered_ipg.to_csv(outfile, header=True, sep="\t", index=False)


    # get the contigs to download
    contigs_to_download = [item[1] for item in parsed_ipg]
    contigs_to_download = list(set([item for sublist in contigs_to_download for item in sublist]))
    # write contigs to download to file
    with open(f"{outdir}/{run_name}/1_ipg/{cl_id}.contigs", "w") as fout:
        #fout.write("\n".join(contigs_to_download))
        for contig in contigs_to_download:
            fout.write(f"{contig}\n")

    return contigs_to_download

def download_contigs(contigs_to_download, already_down_contigs, cl_id, outdir, run_name):
    '''
    '''

    # remove from the list contigs already downloaded
    # use set.difference() for this
    contigs_to_download_set = set(contigs_to_download)
    contigs_to_download = contigs_to_download_set.difference(already_down_contigs)
    contigs_to_download = list(contigs_to_download)
    print(cl_id, len(contigs_to_download_set), len(contigs_to_download))
    print(contigs_to_download)

    outfile = f"{outdir}/{run_name}/2_fasta/{cl_id}.fasta"
    with open(outfile, "w") as fout:
        seqs = str()
        for i in range(0, len(contigs_to_download), 150):
            batch = contigs_to_download[i:i + 150]
            handle = Entrez.efetch(db="nuccore", id=batch, rettype="fasta", retmode="text")
            seqs += handle.read()
            handle.close()

        fout.write(seqs)

    # touch mock file
    os.system(f"touch {outdir}/{run_name}/2_fasta/{cl_id}.done")

def process_mrca_file(already_down_contigs, run_name, mrca_file):
    '''
    '''

    cl_id = os.path.basename(mrca_file).split(".")[0]

    outdir = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial"

    # get Bacteria proteins from polytomies and sister clades
    bacteria_proteins = get_proteins_to_download(mrca_file)
    print(cl_id, len(bacteria_proteins))

    # download proteins info from IPG
    ipg_file = f"{outdir}/{run_name}/1_ipg/{cl_id}.ipg"
    if not os.path.isfile(ipg_file):
        ipg_table = download_ipg(bacteria_proteins, cl_id, outdir, run_name)
    else:
        ipg_table = True

    # if ipg was donwload, it is stored under ipg_table and has been written to file.
    # download genomes then
    if ipg_table and os.path.getsize(ipg_file) != 0:
        # check if the .done mock file has been generated
        if not os.path.isfile(f"{outdir}/{run_name}/2_fasta/{cl_id}.done"):
            contigs_to_download = parse_ipg(cl_id, outdir, run_name)
            download_contigs(contigs_to_download, already_down_contigs, cl_id, outdir, run_name)



def main():

    args = parse_args()


    Entrez.email   = "dcarrillo.bioinf@gmail.com"
    Entrez.api_key = "b6db7fece605d37fcabd4b93749d2e46aa09"


    # list mrca files
    mrca_dir = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/0_mrcas"
    mrca_files = glob.glob(f"{mrca_dir}/*.txt")

    # prepare outdirs
    root_outdir = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial"
    os.makedirs(f"{root_outdir}/{args.run_name}/1_ipg/", exist_ok=True)
    os.makedirs(f"{root_outdir}/{args.run_name}/2_fasta/", exist_ok=True)


    # get which contigs were already downloaded
    already_down_contigs = get_already_downloaded_contigs()


    part = partial(process_mrca_file, already_down_contigs, args.run_name)
    pool = multiprocessing.Pool(processes=13)
    pool.map(part, mrca_files)
    pool.close()
    pool.join()

    # for mrca_file in mrca_files:
    #     process_mrca_file(already_down_contigs, args.run_name, mrca_file)



if __name__ == "__main__":
    main()
