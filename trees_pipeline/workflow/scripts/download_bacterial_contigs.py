from Bio import Entrez
import pandas as pd
from math import ceil
import os

Entrez.email   = snakemake.config['entrez_email']
Entrez.api_key = snakemake.config['entrez_api_key']

# read summary file with the hits
summary_file = f"{snakemake.config['summaries_dir']}/{snakemake.wildcards.cl.replace('_ncbi', '')}.summary"
summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
# keep only Bacteria
bacteria_summary_df = summary_df[(summary_df["included"] == True) & (summary_df["superkingdom"] != "Viruses")]

# download ipg tables in batches of 200
to_download = bacteria_summary_df.index.to_list()
# check if there are Bacteria hits
check = False
if to_download:
	ipg_table = str()
	for i in range(0, len(to_download), 200):
		batch = to_download[i:i + 200]
		handle = Entrez.efetch(db="protein", id=batch, rettype="ipg", retmode="text")
		ipg_table += handle.read().decode("utf-8")
		handle.close()

	with open(snakemake.output.ipg, "w") as fout:
		if ipg_table:
			check = True
			fout.write(ipg_table)

# if there were Bacteria, read the .ipg table and download the contigs
if check:
	lines = [line.strip().split("\t") for line in open(snakemake.output.ipg, "r").readlines()[1:]]
	to_download = list()
	done = list()
	# process RefSeq hits
	for line in lines:
		if line[1] == "RefSeq":
			to_download.append(line[2])
			done.append(line[2].replace("NZ_", ""))
	# process GenBank hits
	for line in lines:
		if line[1] == "INSDC" and line[2] not in done:
			to_download.append(line[2])
			done.append(line[2])


	to_download = list(set(to_download))
	print(snakemake.wildcards.cl, len(to_download))
	with open(snakemake.output.fasta, "w") as fout:
		seqs = str()
		for i in range(0, len(to_download), 150):
			batch = to_download[i:i + 150]
			#print(batch)
			handle = Entrez.efetch(db="nuccore", id=batch, rettype="fasta", retmode="text")
			seqs += handle.read()
			#fout.write(batch_data)
			handle.close()
			#print("batch.gb created")

		fout.write(seqs)

	# write to file all the contigs IDs that should have been downloaded
	with open(snakemake.output.to_download, "w") as fout:
		fout.write("\n".join(to_download))

else:
	os.system(f"touch {snakemake.output.fasta}")
	os.system(f"touch {snakemake.output.ipg}")
	os.system(f"touch {snakemake.output.to_download}")



# # read sister clade with Bacteria hits
# if os.path.getsize(snakemake.input[0]) != 0:
# 	mrcas_sisters_df = pd.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)
# 	mrcas_sisters_df = mrcas_sisters_df.fillna("NA")
# 	# get bacterial contigs to download
# 	if not mrcas_sisters_df.empty:
# 		check = True
# 		polytomies = [polytomy.split(",") for polytomy in mrcas_sisters_df["polytomies"].values.tolist()]
# 		sister_bacteria = [sister.split(",") for sister in mrcas_sisters_df["bacteria_seqs_sister"].values.tolist()]
# 		all = polytomies + sister_bacteria
# 		flat_list = [item for sublist in all for item in sublist]
# 		to_download = list(set(flat_list))
#
# 		if "NA" in to_download:
# 			to_download.remove("NA")
#
# 		print(snakemake.wildcards, to_download)
#
# 		# Download genomes in batches of 200
# 		ipg_table = str()
# 		for i in range(0, len(to_download), 200):
# 			batch = to_download[i:i + 200]
# 			handle = Entrez.efetch(db="protein", id=batch, rettype="ipg", retmode="text")
# 			ipg_table += handle.read().decode("utf-8")
# 			handle.close()
#
# 		with open(snakemake.output.ipg, "w") as fout:
# 			if ipg_table:
# 				check = True
# 				fout.write(ipg_table)
# 			else:
# 				check = False
#
#
# 	else:
# 		check = False
#
#
#
# else:
# 	check = False


# if check:
# 	df = pd.read_csv(snakemake.output.ipg, sep="\t", header=0)
# 	# process the "Nucleotide Accession" table to get unique identifiers to download
# 	to_download = list(set([accession.replace("NZ_", "") for accession in df["Nucleotide Accession"]]))
#
# 	with open(snakemake.output.fasta, "w") as fout:
# 		seqs = str()
# 		for i in range(0, len(to_download), 150):
# 			batch = to_download[i:i + 150]
# 			handle = Entrez.efetch(db="nuccore", id=batch, rettype="fasta", retmode="text")
# 			seqs += handle.read()
# 			#fout.write(batch_data)
# 			handle.close()
# 			#print("batch.gb created")
#
# 		fout.write(seqs)
#
# else:
# 	os.system(f"touch {snakemake.output.fasta}")
# 	os.system(f"touch {snakemake.output.ipg}")















# # keep only Bacteria
# ncbi_summary_df_bacteria = ncbi_summary_df[ncbi_summary_df["superkingdom"] == "Bacteria"]
#
# if not ncbi_summary_df_bacteria.empty:
#     bacterial_hits = ncbi_summary_df_bacteria.index.to_list()
#     ipg_table = str()
#     # split in batches of 200
#     for i in range(0, len(bacterial_hits), 200):
#         batch = bacterial_hits[i:i + 200]
#         handle = Entrez.efetch(db="protein", id=batch, rettype="ipg", retmode="text")
#         ipg_table += handle.read().decode("utf-8")
#         handle.close()
#
#     with open(snakemake.output.ipg, "w") as fout:
#         if ipg_table:
#             check = True
#             fout.write(ipg_table)
#         else:
#             fout.write("No contigs were downloaded.\n")
#
# else:
#     check = False
#     with open(snakemake.output.ipg, "w") as fout:
#         fout.write("No bacterial hits for the CL.\n")
#     with open(snakemake.output.fasta, "w") as fout:
#         fout.write("No bacterial hits for the CL.\n")
#
# if check:
#     df = pd.read_csv(snakemake.output.ipg, sep="\t", header=0)
#     # process the "Nucleotide Accession" table to get unique identifiers to download
#     to_download = list(set([accession.replace("NZ_", "") for accession in df["Nucleotide Accession"]]))
#
#     with open(snakemake.output.fasta, "w") as fout:
#         seqs = str()
#         for i in range(0, len(to_download), 150):
#             batch = to_download[i:i + 150]
#             handle = Entrez.efetch(db="nuccore", id=batch, rettype="fasta", retmode="text")
#             seqs += handle.read()
#             #fout.write(batch_data)
#             handle.close()
#             #print("batch.gb created")
#
#         fout.write(seqs)
