from Bio import SeqIO
import pandas as pd
import numpy as np


def parse_contigs_fasta():
    '''
    Parse all the downloaded contigs and assign a category:
    - ok
    - short, if the contig is shorter than 1.5Kb
    '''

    # parse all contigs, mark them as short if necessary
    all_contigs_records = {record.id:record for record in SeqIO.parse(snakemake.input.all_contigs_fasta, "fasta")}
    contigs_status = {record:"ok" for record in all_contigs_records}
    for record_id, record in all_contigs_records.items():
        if len(record.seq) < 1500:
            contigs_status[record_id] = "short"

    # parse filtered contigs fasta
    filt_contigs_records = {record.id:record for record in SeqIO.parse(snakemake.input.filt_contigs_fasta[0], "fasta")}

    return filt_contigs_records, all_contigs_records, contigs_status

def get_ipgs_for_prots(prots_ids):
    '''
    Parses the ipg table and assign the unique identifiers associate with each protein
    '''

    prots_ipgs = {prot_id:list() for prot_id in prots_ids}

    # iterate the ipg table looking for the unique identifiers for each protein
    ipg_lines = [line.strip().split("\t") for line in open(snakemake.input.ipg, "r").readlines()]
    for line in ipg_lines:
        if line[6] in prots_ipgs:
            prots_ipgs[line[6]] += [line[0]]

    # do the inverse, store in a dict k=ipg_id v=prot_id
    ipgs_prots = dict()
    for prot, ipg_ids in prots_ipgs.items():
        for ipg_id in ipg_ids:
            ipgs_prots[ipg_id] = prot

    return prots_ipgs, ipgs_prots

def get_contigs_ipgs(filt_contigs_records, ipgs_prots):
    '''
    For each ipg, store the contigs associated with it.
    '''

    ipg_lines = [line.strip().split("\t") for line in open(snakemake.input.ipg, "r").readlines()]
    # remove .ipg lines containing contigs that were not downloaded
    ipg_lines = [line for line in ipg_lines if line[2] in filt_contigs_records]
    # iterate the .ipg lines while storing the contigs downloaded for each bacteria protein
    ipgs_contigs = {ipg_id:dict() for ipg_id in ipgs_prots}
    for line in ipg_lines:
        if line[2] in ipgs_contigs[line[0]]:
            ipgs_contigs[line[0]][line[2]] += [[int(line[3]), int(line[4])]]
        else:
            ipgs_contigs[line[0]][line[2]] = [int(line[3]), int(line[4])]

    return ipgs_contigs

def parse_vs2():
    '''

    '''

    # read scores table
    scores_file = snakemake.input.vs2[0].replace("-boundary.tsv", "-score.tsv")
    lines = [line.strip().split("\t") for line in open(scores_file, "r").readlines()[1:]]

    contigs_vs2 = {line[0].split("|")[0]:dict() for line in lines}
    for line in lines:
        contig = line[0].split("|")[0]
        contigs_vs2[contig][line[0]] = [float(line[3])]

    # read boundaries file and store start and stop positions in the dict
    lines = [line.strip().split("\t") for line in open(snakemake.input.vs2[0], "r").readlines()[1:]]
    for line in lines:
        contigs_vs2[line[0]][line[-1]] += [int(line[3]), int(line[4])]

    return contigs_vs2

def check_protein_in_prophage(prot_start, prot_end, prophage_start, prophage_end):
    if prophage_start < prot_start < prophage_end or prophage_start < prot_end < prophage_end:
        return True
    else:
        return False



# read summary file to know which proteins should I create an entry in the table for
summary_file = f"{snakemake.config['summaries_dir']}/{snakemake.wildcards.cl.replace('_ncbi', '')}.summary"
summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
# keep only Bacteria
bacteria_summary_df = summary_df[(summary_df["included"] == True) & (summary_df["superkingdom"] != "Viruses")]
prots_ids = bacteria_summary_df.index.to_list()


# get the relation between prot_ids and ipg_ids
prots_ipgs, ipgs_prots = get_ipgs_for_prots(prots_ids)

# parse downloaded and filtered contigs noting which ones are short
filt_contigs_records, all_contigs_records, contigs_status = parse_contigs_fasta()

# associate contigs to their ipgs
ipgs_contigs = get_contigs_ipgs(all_contigs_records, ipgs_prots)

# parse vs2 results, store in a dict with the prophage regions for each contig
contigs_vs2 = parse_vs2()



# init a dataframe with the bacteria proteins as index
columns = ["ipg_id", "prophage", "discard", "n_contigs", "n_short", "n_prophage", "contigs", "short_contigs", "prophage_contigs"]
final_df = pd.DataFrame(0, index=prots_ids, columns=columns)
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

        # iterate the contigs associated with each ipg
        for contig_id, prot_coords in ipgs_contigs[ipg_id].items():

            # first check it is not a "short" contig
            if contigs_status[contig_id] != "short":
                check = False
                if contig_id in contigs_vs2:
                    for prophage in contigs_vs2[contig_id]:
                        score = contigs_vs2[contig_id][prophage][0]
                        prophage_start = contigs_vs2[contig_id][prophage][1]
                        prophage_end = contigs_vs2[contig_id][prophage][2]

                        if check_protein_in_prophage(prot_coords[0], prot_coords[1], prophage_start, prophage_end):
                            check = True
                    if check:
                        final_df.loc[prot_id, "n_prophage"] += 1
                        final_df.loc[prot_id, "prophage"] = True
                        prophage_contigs.append(contig_id)

        # write prophage contigs, if any
        final_df.loc[prot_id, "prophage_contigs"] = ",".join(prophage_contigs)


    # No contigs were downloaded for the protein
    else:
        final_df.loc[prot_id, "discard"] = True


print(final_df)
