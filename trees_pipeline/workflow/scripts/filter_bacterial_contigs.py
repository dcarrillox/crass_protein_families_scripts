import pandas as pd
from Bio import SeqIO
import os


# read .fasta downloaded contigs to know which ones where downloaded
contigs_records = {record.id:record for record in SeqIO.parse(snakemake.input.fasta, "fasta")}

# read .ipg file and store in a dict k=uniq_id  v=[contig1, contig2...], keep only previously downloaded contigs
lines = [line.strip().split("\t") for line in open(snakemake.input.ipg, "r").readlines()[1:]]
if lines:
    ids_contigs = dict()
    for line in lines:
        if line[2] in contigs_records:
            if line[0] not in ids_contigs:
                ids_contigs[line[0]] = [line[2]]
            else:
                ids_contigs[line[0]] += [line[2]]


    # iterate the ids and their associated contigs. If there are more than 5, keep
    # the top 5 longest
    to_write_ids = list()
    for id, contigs in ids_contigs.items():
        #print(contigs)
        if len(contigs) > 5:
            # get length of the contigs, store in a list of tuples
            contigs_length = [(contig,len(contigs_records[contig].seq)) for contig in contigs]
            sorted_by_length = sorted(contigs_length, key=lambda tuple: int(tuple[1]), reverse=True)
            #print(sorted_by_length)
            to_write_ids += [tuple[0] for tuple in sorted_by_length[:5]]
        else:
            to_write_ids += contigs

    # get the records of the contigs
    to_write_ids = list(set(to_write_ids))
    to_write_records = [contigs_records[contig] for contig in to_write_ids if len(contigs_records[contig].seq) > 1500]

    os.makedirs("results/4_prophages/", exist_ok = True)
    if to_write_records:
        with open(snakemake.output[0], "w") as fout:
            SeqIO.write(to_write_records, fout, "fasta")
    else:
        os.system(f"touch {snakemake.output}")

else:
    os.system(f"touch {snakemake.output}")
