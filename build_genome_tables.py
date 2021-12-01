'''
The input is the directory with all the .faa, one per genome. And the file
"table_OGs_protein_names.txt" from Broccoli step3. Notice that, to show also the
circdup proteins in the tables I don't use the "5_final_faa_circdup_filtered" faa
files but the "3_final_annotation_formatted".
'''

import argparse
import multiprocessing
from functools import partial
import pandas as pd
import numpy as np
from Bio import SeqIO
import glob
import os

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--input_faa_directory',
                               dest='in_dir',
                               required=True,
                               help='directory with final annotation .faa file for each genome.'
                               'I use "2_Prodigal/3_final_annotation_formatted"'
                               )
    requiredArgs.add_argument('-gff', '--input_gff_directory',
                               dest='gff_dir',
                               required=True,
                               help='directory with the GFF prodigal annotation. I use '
                               '"2_Prodigal/3_final_annotation_formatted"'
                               )

    requiredArgs.add_argument('-cls', '--cls_file',
                               dest='cls_file',
                               required=True,
                               help='tabular file with the cl_id for each protein'
                               )
    requiredArgs.add_argument('-scls', '--scls_file',
                               dest='scls_file',
                               required=False,
                               help='output file from mcl containing the SCLs groups'
                               )



    requiredArgs.add_argument('-c', '--coding_file',
                               dest='coding_file',
                               required=True,
                               help='"final_annotation_coding.txt" to get the coding information'
                               )


    requiredArgs.add_argument('-o', '--output_dir',
                               dest='out_dir',
                               required=True,
                               help='directory to put the genome tables, one file per genome'
                               )


    return parser.parse_args()


def CLs_to_dict(cls_file):
    '''
    '''
    lines = [line.strip().split("\t") for line in open(cls_file).readlines()]
    cls_dict = {line[1]:line[0] for line in lines}

    return cls_dict

def SCLs_to_dict(scls_file):
    '''
    For each CL in a SCL, store k=CL  v=SCL
    '''
    lines = [line.strip().split("\t") for line in open(scls_file).readlines()]
    # get how many number should be in the identifier
    n = len(str(len(lines)))

    cont = 1
    scls_dict = dict()
    for line in lines:
        id = str(cont).zfill(n)
        scl_id = f"cl_s_{id}"

        for cl in line:
            scls_dict[cl] = scl_id

        cont += 1

    print(scls_dict)
    return scls_dict


def parse_gff(gff_dir):
    '''
    It reads the gff files in the provided directory and returs a dictionary with
    the start, stop, strand and partial information for each gene.
    k=genome-gene, v=[start, stop, strand, partial]
    '''

    # get the files
    gff_files = glob.glob(f"{gff_dir}/*.gff")

    gff_parsed = dict()
    for gff_file in gff_files:
        genome = os.path.basename(gff_file).split("_prodigal")[0]

        # read the gff file of the genome
        lines = [line.split("\t") for line in open(gff_file).readlines() if not line.startswith("#")]

        for line in lines:
            # get the short_gene_id and create its entry in the dictionary
            gene = line[-1].split("ID=1_")[1].split(";")[0]
            gff_parsed[f"{genome}-{gene}"] = [line[3], line[4], line[6]]

            # if I find an incomplete gene, change the "partial" variable to "True"
            partial = "False"
            if ";partial=10;" in line[-1] or ";partial=01;" in line[-1]:
                partial = "True"
            gff_parsed[f"{genome}-{gene}"].append(partial)

    return gff_parsed

def parse_crAss001_MH675552(gff_dict, nr_annot):
    '''
    Parse the crAss001_MH675552 genes information separately since it does not have
    a gff file. I take this information from the NCBI protein file
    '''
    NCBI_faa_file = "/home/danielc/projects/Bas_phages/2_Prodigal/1_ncbi_references/crAss001_MH675552_NCBI/crAss001_MH675552_NCBI.faa"
    records = SeqIO.parse(NCBI_faa_file, "fasta")

    for record in records:
        start  = str()
        stop   = str()
        strand = str()
        length = len(record.seq)

        short_id = record.description.split("[locus_tag=crAss001_")[1].split("]")[0]
        # I parse the first gene separately since it is the join of two partial genes at the ends
        if short_id == "1":
            tmp = record.description.split("[location=complement(join(102620..102679,")[1].split("))]")[0]
            start  = tmp.split("..")[0]
            stop   =  tmp.split("..")[1]
            strand = "-"
        else:

            if "location=complement(" in record.description:
                tmp = record.description.split("[location=complement(")[1].split(")")[0]
                start  = tmp.split("..")[0]
                stop   =  tmp.split("..")[1]
                strand = "-"
            else:
                tmp = record.description.split("[location=")[1].split("] ")[0]
                start  = tmp.split("..")[0]
                stop   = tmp.split("..")[1]
                strand = "+"

        # add info to gff_dict.
        # notice that there are not partial proteins, so all of them are "False" for that
        gff_dict[f"crAss001_MH675552-{short_id}"] = [start, stop, strand, "False"]

        # add functional info to the NR_dict
        complete_id = f"crAss001_MH675552|{length}|{short_id}"
        protein = record.description.split(" [protein=")[1].split("] ")[0]
        nr_annot[complete_id] = protein

def parse_crAssphage(gff_dict, nr_annot):
    # read the table with start, stop, strand, protein product...
    lines = [line.strip().split("\t") for line in open("/home/danielc/projects/Bas_phages/2_Prodigal/1_ncbi_references/crassphage_NC_024711/GCF_000922395.1_ViralProj259336_feature_table.txt").readlines()[1:]]
    # store in a dict with k=short_id
    product = dict()
    for line in lines:
        if line[0] == "CDS":
            short_id = line[16].split("_gp")[1]
            if short_id.startswith("0"):
                short_id = short_id[1]
            # store start, stop, strand, product
            gff_dict[f"crAssphage-{short_id}"] = [line[7], line[8], line[9], "False"]
            # and also the product for the nr_dict
            product[short_id] = line[13]

    # parse the already formatted .faa file
    records = SeqIO.parse("/home/danielc/projects/Bas_phages/2_Prodigal/1_ncbi_references/crassphage_NC_024711/crAssphage.faa", "fasta")

    for record in records:
        short_id = record.id.split("|")[-1]
        nr_annot[record.id] = product[short_id]

def parse_coding(coding_file):
    '''
    Reads the "final_annotation_coding.txt" with the coding information, per genome,
    and stores this information in a dictionary
    '''
    coding = {line.split("\t")[0]:line.strip().split("\t")[1] for line in open(coding_file).readlines()}

    return coding

def create_table(cls_dict, scls_dict, gff_info, coding_info, yutin_annot, pfam_annot, nr_annot, out_dir, faa_file):
    '''
    Reads the faa file with ALL the proteins of the genome. Then it writes all
    associated information:

    Genome  Coding  Start  End  Strand  Partial?  Yutin  PFAM  NR  Protein_id  CL   SCL
      0       1       2     3      4        5       6      7   8       9       10    11
    '''

    # read faa file
    records = SeqIO.parse(faa_file, "fasta")

    # iterate sequences and add to table
    table = list()
    for record in records:
        # get available information from the protein_id (protein header)
        split = record.id.split("|")
        genome = split[0]
        short  = split[-1]
        start  = gff_info[f"{genome}-{short}"][0]
        end    = gff_info[f"{genome}-{short}"][1]
        strand = gff_info[f"{genome}-{short}"][2]
        partialness = gff_info[f"{genome}-{short}"][3]
        coding_genome = coding_info[genome]

        # check for yutin
        full_length = False
        if record.id in yutin_annot:
            for annot in yutin_annot[record.id]:
                #print(annot)
                if annot[3] == "full_length":
                    yutin = annot[2]
                    full_length = True
                    break
            # if there are only domain hits, pick the first hit
            if not full_length:
                yutin = yutin_annot[record.id][0][2]
        else:
            yutin = ""

        # check for pfam
        if record.id in pfam_annot:
            pfam = pfam_annot[record.id]
        else:
            pfam = ""

        # check for nr
        if record.id in nr_annot:
            nr= nr_annot[record.id]
        else:
            nr = ""

        # check for OG_id
        if record.id in cls_dict:
            cl_id = cls_dict[record.id]
        else:
            cl_id = ""

        # check for SCL
        if cl_id in scls_dict:
            scl_id = scls_dict[cl_id]
        else:
            scl_id = ""

        # merge all the information and add to table
        tow = [genome, coding_genome, start, end, strand, partialness, yutin, pfam, nr, record.id, cl_id, scl_id] # removed circdup
        table.append(tow)


    outfile = f"{out_dir}/" + os.path.basename(faa_file).replace(".faa", ".table")
    with open(outfile, "w") as fout:
        header = ["genome", "coding", "start", "end", "strand", "partial", "yutin", "pfam", "nr", "protein_id","CL","SCL"] # removed circdup
        fout.write("\t".join(header) + "\n")

        for line in table:
            fout.write("\t".join(line) + "\n")



def main():

    args = parse_args()

    ## To generate the 0_raw tables, uncomment this block to provide empty dicts for functional
    ## annotation. To generate the 1_function, run the block below below where the parsed
    ## annotation results are provided
    yutin_annot = dict()
    pfam_annot = dict()
    nr_annot = dict()

    lines = [line.strip().split("\t") for line in open("/home/danielc/projects/Bas_phages/2_Function/parsed_yutin_all.txt").readlines()]
    # init the dict with the keys (protein_ids)
    yutin_annot = {line[0]:list() for line in lines}
    for line in lines:
        yutin_annot[line[0]] += [line]
    #yutin_annot = {line.split("\t")[0]:line.strip().split("\t") for line in open("/home/danielc/projects/Bas_phages/2_Function/parsed_yutin_all.txt").readlines()}
    #print(yutin_annot)
    pfam_annot = {line.split("\t")[0]:line.strip().split("\t")[1] for line in open("/home/danielc/projects/Bas_phages/2_Function/parsed_pfam.txt").readlines()}
    # nr_annot = dict()



    # read broccoli table relating each protein to an OG
    CL_annot = CLs_to_dict(args.cls_file)
    if args.scls_file:
        SCL_annot = SCLs_to_dict(args.scls_file)
    else:
        SCL_annot = dict()

    # parse gff to get start, stop, strand and partial info
    gff_info = parse_gff(args.gff_dir)
    #print(gff_info)
    # notice that crAss001_MH675552 has not gff, so I have to read the info from the
    # NCBI protein file and add it to the gff_dict
    # add also the protein info to the nr_annot dict
    parse_crAss001_MH675552(gff_info, nr_annot)
    parse_crAssphage(gff_info, nr_annot)

    # read the
    coding_info = parse_coding(args.coding_file)

    # Now get faa files. I use these because they contain all the proteins, including
    # those that were not assigned to an OG.
    faa_files = glob.glob(args.in_dir + "/*.faa")
    # remove FUFK010039141 genome
    # faa_files.remove(f"{args.in_dir}/FUFK010039141.faa")
    # print(len(faa_files))

    # create partial function
    part = partial(create_table, CL_annot, SCL_annot, gff_info, coding_info, yutin_annot, pfam_annot, nr_annot, args.out_dir) # circdup_proteins removed

    pool = multiprocessing.Pool(processes=18)
    shared_content = pool.map(part, faa_files)
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
