'''
This script parses the hmmsearch results of searching proteins with profiles. It
returns two files:
- All the results (parsed-all)
- Protein hit to its CL (parsed-CLhit)
'''

import argparse
import glob
import os
from Bio import SeqIO, SearchIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_file',
                               dest='hmmsearch_file',
                               required=True,
                               help='hmmer3-text file with the hmmsearch results'
                               )

    requiredArgs.add_argument('-f', '--faa_dir',
                               dest='faa_dir',
                               required=True,
                               help='directory with the CL.faa files. I parse them '
                               'to know to which CL the proteins were assigned'
                               )
    requiredArgs.add_argument('-m', '--multi_faa',
                               dest='multi_faa',
                               required=True,
                               help='faa file with the multi proteins'
                               )

    requiredArgs.add_argument('-s', '--summary_file',
                               dest='summary_file',
                               required=True,
                               help='file "separate_cls.txt" with the info of every CL '
                               'being single or multi'
                               )

    requiredArgs.add_argument('-o', '--out_file_prefix',
                               dest='out_file_prefix',
                               required=True,
                               help='prefix of the two outfiles.'
                               )

    return parser.parse_args()


def parse_hmmsearch(hmmsearch_file):
    # read hmmer3-text file
    records = SearchIO.parse(hmmsearch_file, "hmmer3-text")

    # store hits in two dicts:
    # CLs_hits: k=CL, v=[[protein hit 1], [protein hit 2] ...]
    # proteins_hits: k=protein_id, v=hit_CL
    CLs_hits = dict()
    proteins_hits = dict()
    for record in records:
        CL_id = record.id
        CLs_hits[CL_id] = list()
        for hit in record.hits:
            if hit.is_included:
                if hit.id in proteins_hits:
                    proteins_hits[hit.id].append(CL_id)
                else:
                    proteins_hits[hit.id] = [CL_id]

                for hsp in hit.hsps:
                    if hsp.is_included:
                        CLs_hits[CL_id].append([hit.id,
                                                 str(hsp.evalue),
                                                 str(hsp.bitscore),
                                                 str(record.seq_len),
                                                 str(hsp.env_start),
                                                 str(hsp.env_end),
                                                 str(hsp.query_start),
                                                 str(hsp.query_end)
                                                 ])

    return CLs_hits, proteins_hits

def get_multi_dicts(multi_faa, faa_dir, summary_file):

    # store multi proteins ids in a set
    multi_proteins = {record.id for record in SeqIO.parse(multi_faa, "fasta")}

    # read "separate_cls.txt" to know which CLs are multi.
    multi_CLs = [line.split("\t")[0] for line in open(summary_file).readlines()[1:] if line.strip().split("\t")[1] == "multi"]

    proteins_CLs = dict()
    CLs_proteins = {CL:list() for CL in multi_CLs}

    faa_files = glob.glob(f"{faa_dir}/*.faa")
    for faa_file in faa_files:
        CL_id = os.path.basename(faa_file).split(".faa")[0]
        if CL_id in multi_CLs:
            with open(faa_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in multi_proteins:
                        proteins_CLs[record.id] = CL_id
                        CLs_proteins[CL_id] += [record.id]

    return proteins_CLs, CLs_proteins


def main():

    args = parse_args()

    print("Parsing hmmsearch file...")
    CLs_hits, proteins_hits = parse_hmmsearch(args.hmmsearch_file)
    #print(CLs_hits)

    print("Writing the 'all' version...")
    with open(f"{args.out_file_prefix}.parsed-all", "w") as fout:
        header = ["protein", "evalue", "bitscore", "qlen", "qstart", "qend", "hstart", "hend"]
        fout.write("\t".join(header) + "\n")
        # I want to sort the hits by genome, and then by n_gene
        for CL, proteins in CLs_hits.items():
            # first get the uniq genomes per CL
            genomes = list(set([protein[0].split("|")[0] for protein in proteins]))
            genomes.sort()
            # then store the hits per genome and sort them by n_gene
            for genome in genomes:
                tow = list()
                for protein in proteins:
                    if protein[0].split("|")[0] == genome:
                        tow.append(protein)

                ordered = sorted(tow, key=lambda protein: int(protein[0].split("|")[-1]))

                for protein in ordered:
                    fout.write(f"{CL}\t" + "\t".join(protein) + "\n")


    print("Parsing multi CLs faa files...")
    proteins_CLs, CLs_proteins = get_multi_dicts(args.multi_faa, args.faa_dir, args.summary_file)


    print("Writing the 'CL' version...")
    # iterate the hmmsearch results
    with open(f"{args.out_file_prefix}.parsed-CLhit", "w") as fout:
        header = ["protein", "evalue", "bitscore", "qlen", "qstart", "qend", "hstart", "hend"]
        fout.write("\t".join(header) + "\n")

        # iterate ALL the CLs that were found to be duplicated (I read them from
        # the summary file and call them "ref") and identify those that didn't get
        # any hit

        # If the CL didn't get any hit...
        CLs_not_hits =list()
        for CL_ref, proteins_ref in CLs_proteins.items():
            if f"{CL_ref}" not in CLs_hits:
                CLs_not_hits.append(CL_ref)
        # write those CLs and their protein to file
        for CL_no_hit in CLs_not_hits:
            for protein in CLs_proteins[CL_no_hit]:
                fout.write(f"{CL_no_hit}\t{protein}\tNo hit\n")

            # write empty line, useful when manually reviewing results
            fout.write("\n")

        # Then iterate the CLs that got hits
        for CL, proteins in CLs_hits.items():
            all_proteins = list()
            # check that the hit CL matches with the protein CL
            matched = [protein for protein in proteins if f"{proteins_CLs[protein[0]]}" == CL]
            matched_ids = [protein[0] for protein in proteins if f"{proteins_CLs[protein[0]]}" == CL]
            print(CL, matched_ids)
            # get those proteins whose CL_hit didn't match
            unmatched = list()
            for protein_id, CL_id in proteins_CLs.items():
                if f"{CL_id}" == CL:
                    if protein_id not in matched_ids:
                        unmatched.append(protein_id)

            for protein_id in unmatched:
                if protein_id in proteins_hits:
                    all_proteins.append([protein_id, "Other hits", ",".join(proteins_hits[protein_id])])
                else:
                    all_proteins.append([protein_id, "No hit"])

            # merge matched and unmatched proteins
            all_proteins += matched

            # I want to sort the hits by genome, and then by n_gene
            # first get the uniq genomes per CL
            genomes = list(set([protein[0].split("|")[0] for protein in all_proteins]))
            genomes.sort()
            # then store the hits per genome and sort them by n_gene
            for genome in genomes:
                tow = list()
                for protein in all_proteins:
                    if protein[0].split("|")[0] == genome:
                        tow.append(protein)

                ordered = sorted(tow, key=lambda protein: int(protein[0].split("|")[-1]))

                for protein in ordered:
                    fout.write(f"{CL}\t" + "\t".join(protein) + "\n")

            # write empty line, useful when manually reviewing results
            fout.write("\n")

if __name__ == "__main__":
    main()
