import glob, os, argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hhr_dir',
                               dest='hhr_dir',
                               required=True,
                               help='directory with the .hhsearch files'
                               )
    requiredArgs.add_argument('-m', '--msa_dir',
                               dest='msa_dir',
                               required=True,
                               help='directory with the .mafft files, used '
                               'to make get the relation sequence-PF_id.'
                               )

    # requiredArgs.add_argument('-omcl', '--out_mcl_file',
    #                            dest='out_mcl_file',
    #                            required=False,
    #                            help='file with the table already formatted for running '
    #                            'MCL. Three columns: CL-1, CL-2, prob*coverage'
    #                            )


    return parser.parse_args()


def fix_hit_line(line):
    # split by space and remove empty fields
    fields = line.split(" ")
    while "" in fields:
        fields.remove("")
    # check if last field has the two last columns mixed
    # if so, split it in two columns and remove the parenthesis...
    if len(fields) != 11:
        mixed = fields[-1].replace(")", "").split("(")
        fields[-1] = mixed[0]
        fields.append(mixed[1])
    # ...if not, just replace the parenthesis in the last columnn
    else:
        fields[-1] = fields[-1].replace("(", "").replace(")", "")

    return fields

def hits_table_to_dict(hits_lines, query_length):
    No_hits = dict()
    # query self alignment length
    query_self_length = fix_hit_line(hits_lines[0])[-1]

    # hits lines
    for line in hits_lines:
        fields = fix_hit_line(line)

        no_hit = fields[0]
        evalue = fields[3]
        prob   = fields[2]
        pvalue = fields[4]
        score  = fields[5]
        cols   = fields[7]
        qstart = fields[8].split("-")[0]
        qend   = fields[8].split("-")[1]
        tstart = fields[9].split("-")[0]
        tend   = fields[9].split("-")[1]
        qlen   = query_length
        #qselfl = "{:.2f}".format(int(query_self_length)/int(query_length))

        No_hits[no_hit] = [evalue, prob, pvalue, score, cols,
                           qstart, qend, tstart, tend] #, qlen, qselfl]

    return No_hits

def seqs_to_cls_ids(align_lines, seqs_cls_ids):
    No_seqs = dict()
    for i, line in enumerate(align_lines):
        if line.startswith("No ") and align_lines[i+1].startswith(">"):
            No_hit = line.split("No ")[1].strip()
            seq_id = align_lines[i+1].replace(">", "").strip()
            No_seqs[No_hit] = seqs_cls_ids[seq_id]

    return No_seqs

def tab_2_mcl(tab_dir, out_mcl_file):
    # list .tab files
    tab_files = glob.glob(f'{tab_dir}/*.tab')

    pfs_evalues = dict()
    # store the evalues for each PF pair in a dict with k=pf_id  v={pf_1:evalue, pf_2:evalue}
    for tab_file in tab_files:
        # get the id of the pf from the file's name
        pf_id = os.path.basename(tab_file).replace(".tab", "")
        # read the file
        lines = [line.strip().split('\t') for line in open(tab_file).readlines()[1:]] # discard header
        # create the dictionary, evalue cutoff of 0.01
        # NB that there might be several hits from the same PF to another. In that case,
        # I take the longest hit. So, store the cols column too. Also, don't account for self hits
        pfs_evalues[pf_id] = {line[1]:list() for line in lines if float(line[2]) < 0.01 and line[0] != line[1]}
        for line in lines:
            if float(line[2]) < 0.01 and line[0] != line[1]:
                # check if there is a previous hit
                # if so, compare number of aligned cols
                if pfs_evalues[pf_id][line[1]]:
                    if int(line[6]) > pfs_evalues[pf_id][line[1]][1]:
                        pfs_evalues[pf_id][line[1]] = [float(line[2]), int(line[6])]
                # if there is no previous hit
                else:
                    pfs_evalues[pf_id][line[1]] = [float(line[2]), int(line[6])]

    # iterate the dict
    pfs_evalues_final = list()
    for pf_query, pfs_targets in pfs_evalues.items():
        #print(pf_query)
        for pf_target, evalue_cols in pfs_targets.items():
            # check that the hit was reciprocal and that it has not been processed yet
            if pf_query in pfs_evalues[pf_target] and evalue_cols[0] != "done":
                mean = np.mean([evalue_cols[0], pfs_evalues[pf_target][pf_query][0]])
                pfs_evalues_final.append([pf_query, pf_target, str(mean)])
                pfs_evalues[pf_target][pf_query][0] = "done"

    with open(out_mcl_file, "w") as fout:
        for line in pfs_evalues_final:
            fout.write("\t".join(line) + "\n")



def main():

    args = parse_args()

    # start by reading the MSAs and getting the first sequence, this is the one
    # that hhsearch shows. Associate it with the cl_id
    msas_files = glob.glob(f'{args.msa_dir}/*.mafft')

    seqs_cls_ids = dict()
    for msa_file in msas_files:
        cl_id = os.path.basename(msa_file).replace('.mafft', '')
        # get the first line of the file. It should be the header of the first seq
        lines = [line.strip() for line in open(msa_file).readlines()]
        if lines:
            seq_id = lines[0].replace(">", "")
            seqs_cls_ids[seq_id] = cl_id
        else:
            print(f"{cl_id} looks empty!")


    # get the length of every profile by looking at the second row of the hhr file
    cls_lengths = {os.path.basename(hhr_file).replace(".hhr", ""):int(open(hhr_file).readlines()[1].split("Match_columns")[1].strip())
                    for hhr_file in glob.glob(f"{args.hhr_dir}/*.hhr")}


    # now iterate the hhr files
    hhr_files = glob.glob(f"{args.hhr_dir}/*.hhr")
    for hhr_file in hhr_files:
        # get id of the cl
        query_cl_id = os.path.basename(hhr_file).split(".hhr")[0]

        # get the lines starting from the header of the table
        lines = [line.strip() for line in open(hhr_file).readlines()[8:]]

        # identify the number of hits in the table. For this, look at the last
        # hit at the hits summary on the top part of the hhr file, and get it n identifier (value before hit's name)
        for i, line in enumerate(lines):
            if line == "No 1":
                n_hits = int(lines[i-2].split(" ")[0])

        # save the hits in a dict, k=No_hit  v=[line.split(" ")]
        hits_lines = lines[1:n_hits+1]
        No_hits = hits_table_to_dict(hits_lines, cls_lengths[query_cl_id])

        # parse the alignments section to get the cl_id of the hits
        align_lines = lines[n_hits+2:]
        No_seqs = seqs_to_cls_ids(align_lines, seqs_cls_ids)

        #
        out_file = hhr_file.replace(".hhr", ".tab")
        with open(out_file, "w") as fout:
            header = ["query", "target", "evalue", "prob", "pvalue", "score", "cols", "qstart",
                      "qend", "tstart", "tend", "qlen", "tlen", "qcov", "tcov", "qscore", "tscore"]
            fout.write("\t".join(header) + "\n")

            # write the hits sorted by No (number)
            # check also if the first hit is the cl itself
            for i in range(1, n_hits+1):
                target_cl_id = No_seqs[str(i)]

                # get alingment length in the query and target
                qaln = abs(int(No_hits[str(i)][6]) - int(No_hits[str(i)][5]))
                taln = abs(int(No_hits[str(i)][8]) - int(No_hits[str(i)][7]))

                qcov = "{:.2f}".format(qaln/cls_lengths[query_cl_id])
                tcov = "{:.2f}".format(taln/cls_lengths[target_cl_id])
                #print(float(No_hits[str(i)][1]), qaln, qcov)
                qscore = "{:.2f}".format(float(No_hits[str(i)][1]) * float(qcov))
                tscore = "{:.2f}".format(float(No_hits[str(i)][1]) * float(tcov))

                to_write = [query_cl_id, target_cl_id] + No_hits[str(i)] + [str(cls_lengths[query_cl_id]), str(cls_lengths[target_cl_id]), qcov, tcov, qscore, tscore]
                fout.write("\t".join(to_write) + "\n")


    # create the file formated for MCL clustering
    #tab_2_mcl(args.hhr_dir, args.out_mcl_file)

if __name__ == "__main__":
    main()
