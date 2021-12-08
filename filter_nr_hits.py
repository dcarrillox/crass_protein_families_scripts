'''
For each of the hits in the hmmsearch file, calculates coverages and put it in
a table. Then open the fasta file and get the last common ancestor of the hits.
'''

from Bio import SeqIO, SearchIO
import os, glob, argparse, multiprocessing, time, itertools
from ete3 import NCBITaxa
from functools import partial
import pandas as pd
import numpy as np



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_dir',
                               dest='hmmsearch_dir',
                               required=True,
                               help='directory with the hmmsearch results'
                               )
    requiredArgs.add_argument('-f', '--faa_dir',
                               dest='faa_dir',
                               required=True,
                               help='directory with the .faa files obtained with '
                               '"extract_seqs_NR.py"'
                               )
    requiredArgs.add_argument('-s', '--summary_dir',
                               dest='summary_dir',
                               required=True,
                               help='directory to save the .summary files'
                               )
    requiredArgs.add_argument('-o', '--outdir',
                               dest='outdir',
                               required=True,
                               help='directory to save the .faa files'
                               )

    return parser.parse_args()



def calculate_coverages(hit, cl_len):
    p_len = hit.seq_len
    hit_evalue =  hit.evalue
    cl_starts, cl_ends = list(), list()
    p_starts, p_ends = list(), list()
    for hsp in hit.hsps:
        if hsp.evalue < 0.001:
            cl_starts.append(hsp.query_start)
            cl_ends.append(hsp.query_end)
            p_starts.append(hsp.env_start)
            p_ends.append(hsp.env_end)

    # protein coverage
    if p_starts:
        intervals = [[s,e] for s,e in zip(p_starts, p_ends)]
        intervals.sort(key=lambda interval: interval[0])
        merged = [intervals[0]]
        for current in intervals:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)

        # calculate coverage
        covered_length = float(0)
        for interval in merged:
            covered_length += 1 + (interval[1] - interval[0])

        p_cov = round(float(covered_length/p_len),3)
    else:
        p_cov = 0


    # CL coverage
    if p_starts:
        intervals = [[s,e] for s,e in zip(cl_starts, cl_ends)]
        intervals.sort(key=lambda interval: interval[0])
        merged = [intervals[0]]
        for current in intervals:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)

        # calculate coverage
        covered_length = float(0)
        for interval in merged:
            covered_length += 1 + (interval[1] - interval[0])

        cl_cov = round(float(covered_length/cl_len),3)
    else:
        cl_cov = 0

    return p_len, p_cov, cl_cov

def parse_hmmsearch_file(hmmsearch_dir, cl_id):
    file = f"{hmmsearch_dir}/{cl_id}.domtblout"

    # check size
    mb = os.path.getsize(file)/(1024*1024)
    exclude = ["cl_2650", "cl_3314", "cl_s_526", "cl_s_321"]
    if cl_id not in exclude:

        # process the file
        to_write = list()
        with open(file, "r") as handle:
            for record in SearchIO.parse(file, "hmmsearch3-domtab"):
                if record.id != cl_id:
                    print("wut?")
                # iterate the hits, keep only those with evalue < 0.1
                cl_len = record.seq_len
                for hit in record.hits:
                    if hit.evalue < 0.1:
                        p_len, p_cov, cl_cov = calculate_coverages(hit, cl_len)
                        if p_cov >= 0.4 and cl_cov >= 0.4:
                            include = True
                        else:
                            include = False

                        to_add = [hit.id, include, hit.evalue, hit.bitscore, p_len, cl_len, p_cov, cl_cov]
                        to_write.append(to_add)

        table = pd.DataFrame(to_write, columns=["protein_id", "included", "evalue", "bitscore", "p_len", "cl_len", "p_cov", "cl_cov"])
        table = table.set_index("protein_id")
        return table

    # file too big, skip it and return None
    else:
        print(f"Skipping file {file}, CL is in excluded, {round(mb, 3)} Mb.")
        return None

def all_equal(iterable):
    # https://docs.python.org/3/library/itertools.html#itertools-recipes
    """
    Returns True if all the elements are equal to each other
    """
    g = itertools.groupby(iterable)
    return next(g, True) and not next(g, False)

def get_lcas(table, faa_dir, outdir, cl_id):
    ncbi = NCBITaxa()

    table["lca"] = np.nan
    table["rank"] = np.nan
    table["nprots"] = np.nan

    # insert column for virus annotated as crass-like in the NCBI
    table.insert(1, "ncbi_crass", False)

    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus"]
    table["superkingdom"] = np.nan
    table["phylum"] = np.nan
    table["class"] = np.nan
    table["order"] = np.nan
    table["family"] = np.nan
    table["genus"] = np.nan

    # list to store eukaryotic record.ids
    euk_arc_records_ids = list()

    # count included sequences
    n_included = table.included.sum()
    #print(n_included)
    if n_included > 0:
        # filter table df to keep only included
        included = table.loc[table.included, :]

        # parse faa file, keep only included
        # check that .faa file is present
        faa_file = f"{faa_dir}/{cl_id}.faa"
        if os.path.isfile(faa_file):
            records = [record for record in SeqIO.parse(faa_file, "fasta") if record.id in included.index]
            #print(records)
            #
            for record in records:
                #print(record.id)
                names  = [name.split(" [")[-1][:-1] for name in record.description.split("\x01")]
                #print(names)

                taxids = list()
                for name in names:
                    try:
                        taxid = ncbi.get_name_translator([name])[name][0]
                        taxids.append(taxid)
                    except:
                        pass

                if taxids:
                    lineages = list()
                    for taxid in taxids:
                        try:
                            full_taxonomy = ncbi.get_lineage(taxid)
                            lineages.append(full_taxonomy)
                        except:
                            pass

                    per_level = [i for i in itertools.zip_longest(*lineages)]
                    lca = per_level[0]
                    for i in range(len(per_level)):
                        if all_equal(per_level[i]):
                            lca = per_level[i][0]
                    # Following the example LCA = 3
                    name_dic = ncbi.get_taxid_translator([lca]) # Translate to name
                    name = name_dic[lca] # Get the name
                    rank_dic = ncbi.get_rank([lca]) # Translate to rank
                    rank = rank_dic[lca] # Get the rank

                    table.loc[record.id, "lca"] = name
                    table.loc[record.id, "rank"] = rank
                    table.loc[record.id, "nprots"] = len(taxids)

                    # include the rest of the ranks
                    lca_lineage = ncbi.get_lineage(lca)
                    ranks = ncbi.get_rank(lca_lineage)
                    for taxid_lineage, rank in ranks.items():
                        if rank in target_ranks:
                            table.loc[record.id, rank] = ncbi.get_taxid_translator([taxid_lineage])[taxid_lineage]

                    # check if the crass taxid is in the lineage
                    if 1978007 in lca_lineage:
                        table.loc[record.id, "ncbi_crass"] = True

                    # check it is not eukaryotic or archaeal
                    if table.loc[record.id, "superkingdom"] in ["Eukaryota", "Archaea"]:
                        euk_arc_records_ids.append(record.id)
                        table.loc[record.id, "included"] = False

            # filter out eukaryotic and archaeal
            records_no_euk_arc = [record for record in records if record.id not in euk_arc_records_ids]

            # write .faa file
            out_faa = f"{outdir}/{cl_id}.faa"
            for record in records_no_euk_arc:
                record.description = ""
            with open(out_faa, "w") as fout:
                SeqIO.write(records_no_euk_arc, out_faa, "fasta")

        else:
            print(f"{faa_dir}/{cl_id}.faa was not found. Skipping taxa assessment")

    return table

def process_CL(hmmsearch_dir, faa_dir, outdir, summary_dir, cl_id):
    hmm_stats  = parse_hmmsearch_file(hmmsearch_dir, cl_id)

    # check if the domtblout file was processed by looking if the hmm_stats variable
    # is None. This is what I return if the file was not processed
    if not isinstance(hmm_stats, type(None)):
        taxa_annot = get_lcas(hmm_stats, faa_dir, outdir, cl_id)

        outfile = f"{summary_dir}/{cl_id}.summary"
        taxa_annot.to_csv(outfile, index=True, header=True, sep="\t")




def main():

    args = parse_args()

    # Update the NCBI database if necessary
    #ncbi = NCBITaxa()
    #ncbi.update_taxonomy_database()

    # get done summaries
    done = [os.path.basename(file).split(".")[0] for file in glob.glob(f"{args.summary_dir}/cl_*.summary")]
    print(f"{len(done)} CLs already have a summary file")

    cl_ids = [os.path.basename(file).split(".")[0] for file in glob.glob(f"{args.hmmsearch_dir}/*.domtblout")
              if os.path.basename(file).split(".")[0] not in done]
    print(f"{len(cl_ids)} CLs will be processed")

    time.sleep(3)

    part = partial(process_CL, args.hmmsearch_dir, args.faa_dir, args.outdir, args.summary_dir)
    pool = multiprocessing.Pool(processes=70)
    pool.map(part, cl_ids)
    pool.close()
    pool.join()



if __name__ == "__main__":
    main()
