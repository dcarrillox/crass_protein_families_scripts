from Bio import Entrez
import pandas as pd
from math import ceil

Entrez.email   = snakemake.config['entrez_email']
Entrez.api_key = snakemake.config['entrez_api_key']

# read summary file of NCBI hits
# notice I don't have the _ncbi flag in the summary files, but it is in my wildcard cl
summary_file = f"{snakemake.config['summaries_dir']}/{snakemake.wildcards.cl.replace('_ncbi','')}.summary"
ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
# keep only Bacteria
ncbi_summary_df_bacteria = ncbi_summary_df[ncbi_summary_df["superkingdom"] == "Bacteria"]

if not ncbi_summary_df_bacteria.empty:
    bacterial_hits = ncbi_summary_df_bacteria.index.to_list()
    ipg_table = str()
    # split in batches of 200
    for i in range(0, len(bacterial_hits), 200):
        batch = bacterial_hits[i:i + 200]
        handle = Entrez.efetch(db="protein", id=batch, rettype="ipg", retmode="text")
        ipg_table += handle.read().decode("utf-8")
        handle.close()

    with open(snakemake.output.ipg, "w") as fout:
        if ipg_table:
            check = True
            fout.write(ipg_table)
        else:
            fout.write("No contigs were downloaded.\n")

else:
    check = False
    with open(snakemake.output.ipg, "w") as fout:
        fout.write("No bacterial hits for the CL.\n")
    with open(snakemake.output.fasta, "w") as fout:
        fout.write("No bacterial hits for the CL.\n")

if check:
    df = pd.read_csv(snakemake.output.ipg, sep="\t", header=0)
    # process the "Nucleotide Accession" table to get unique identifiers to download
    to_download = list(set([accession.replace("NZ_", "") for accession in df["Nucleotide Accession"]]))

    with open(snakemake.output.fasta, "w") as fout:
        seqs = str()
        for i in range(0, len(to_download), 150):
            batch = to_download[i:i + 150]
            handle = Entrez.efetch(db="nuccore", id=batch, rettype="fasta", retmode="text")
            seqs += handle.read()
            #fout.write(batch_data)
            handle.close()
            #print("batch.gb created")

        fout.write(seqs)
