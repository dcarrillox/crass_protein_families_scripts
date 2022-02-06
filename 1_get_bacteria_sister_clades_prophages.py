'''

'''

import os, glob, argparse, multiprocessing, time, itertools
from functools import partial
import pandas as pd
import numpy as np
from ete3 import Tree, RectFace


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-n', '--run_name',
                               dest='run_name',
                               required=True,
                               help='name of the analysis (round1, round2...)'
                               )

    return parser.parse_args()



def get_crassvirales_taxonomy_df(taxa_file):
    '''
    Parses the taxonomy and returns a dataframe, index=genomes, 3 columns (family, subfamily, genus)
    '''
    #"/home/danielc/software/github/old_june/tree_functions/files/crassphages_taxonomy_terL_and_new.txt"
    df = pd.read_csv(taxa_file, sep="\t", index_col=0, names=["family","subfamily","genus"])
    df['superkingdom'] = np.where(df['family'].isin(["new", "crass_env"]), "new", "crassvirales")
    return df

def create_target_taxa_groups():
    taxa_groups_raw = ["crassvirales","NA", "new", "Viruses", "bact_crass_polytomy", "prophage"]
    taxa_groups = list()
    # create the target taxa_groups to look for in the tree
    for i in range(2, len(taxa_groups_raw)+1):
        combinations = itertools.combinations(taxa_groups_raw, i)
        # iterate the combinations and retain those with "crassvirales" in them
        for combination in combinations:
            if "crassvirales" in combination:
                taxa_groups.append(list(combination))

    return taxa_groups

def parse_cl_vs2_results(cl_id):
    '''
    '''
    # check if there are vs2 results for this cl_id
    cl_vs2_files = glob.glob(f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/*/3_vs2/2_classified/{cl_id}.prophages")
    cl_vs2_prophages = dict()
    if cl_vs2_files:
        cl_vs2_lines = list()
        for file in cl_vs2_files:
            lines = [line.strip().split("\t") for line in open(file).readlines()[1:]]
            cl_vs2_lines += lines

        # store in a dict if prophage, short or not_down
        for line in cl_vs2_lines:
            # check n_contigs to know if the contig's protein was analyzed
            if line[4] != "0":
                protein_id = line[0]
                # it is a prophage
                if line[2] == "True":
                    cl_vs2_prophages[protein_id] = "prophage"
                if line[4] == line[6]:
                    cl_vs2_prophages[protein_id] = "short"
                if line[4] == line[7]:
                    cl_vs2_prophages[protein_id] = "not_down"

    return cl_vs2_prophages

def resolve_crass_bacteria_polytomies(t):
    '''
    Looks for clades with Crassvirales and Bacteria, all of them with 0 dist
    '''
    checked_leaves = list()
    crass_polytom_leaves = t.search_nodes(dist=0.0, superkingdom="crassvirales")
    for leaf in crass_polytom_leaves:
        if leaf.name not in checked_leaves:
            checked_leaves.append(leaf.name)
            # go one node up, check superkingdom of descendants
            node_up = leaf.up
            for node in node_up.iter_descendants():
                if node.superkingdom == "Bacteria":
                    node.superkingdom = "bact_crass_polytomy"
                    print(f"\t- polytomy bacteria-crass, {node.name}")

    return t

def read_and_annotate_tree(tree_file, ncbi_summary_df, crassvirales_taxa_df, cl_vs2_prophages):
    '''
    Reads the tree of a CL with NCBI and Crassvirales sequences and assigns
    source and taxonomy to them
    '''

    # read tree
    t = Tree(tree_file, format=1)

    # iterate the leaves, assi
    for leaf in t.iter_leaves():
        # check if it is an NCBI hit:
        if leaf.name in ncbi_summary_df.index:
            leaf.add_features(source="ncbi",
                              crass=str(ncbi_summary_df.loc[leaf.name, "ncbi_crass"]),
                              superkingdom=ncbi_summary_df.loc[leaf.name, "superkingdom"],
                              phylum=ncbi_summary_df.loc[leaf.name, "phylum"],
                              #class=ncbi_summary_df.loc[leaf.name, "class"],
                              order=ncbi_summary_df.loc[leaf.name, "order"],
                              family=ncbi_summary_df.loc[leaf.name, "family"],
                              genus=ncbi_summary_df.loc[leaf.name, "genus"]
                              )

            # check if the protein was called as prophage by vs2
            if leaf.name in cl_vs2_prophages:
                leaf.superkingdom = cl_vs2_prophages[leaf.name]
                print(leaf.name, leaf.superkingdom)



        # if not, it is a Crassvirales sequence
        else:
            genome = leaf.name.split("|")[0]
            leaf.add_features(source="crass",
                              crass=str(True),
                              superkingdom=crassvirales_taxa_df.loc[genome, "superkingdom"],
                              # I put the Crassvirales family here
                              phylum=crassvirales_taxa_df.loc[genome, "family"],
                              #class="crass",
                              order="crass",
                              family=crassvirales_taxa_df.loc[genome, "family"],
                              genus=crassvirales_taxa_df.loc[genome, "genus"]
                              )

    return t

def identify_crassvirales_mrcas(t, taxa_groups, start_cont, faces=False):
    '''

    '''

    monophyletic_clades_leaves = dict()
    monophyletic_clades = dict()
    mrcas_ids_to_clean = list()
    cont = 0
    for taxa_group in taxa_groups:
        mrcas = t.get_monophyletic(target_attr="superkingdom", values=taxa_group)
        if mrcas:
            for mrca in mrcas:
                cont += 1
                mrca_id = f"mrca_{cont}"
                # get names of the leaves
                mrca_leaves_names = {leaf.name for leaf in mrca.iter_leaves()}
                if monophyletic_clades_leaves:
                    for target_mrca, target_leaves_names in monophyletic_clades_leaves.items():
                        tmp = target_leaves_names.union(mrca_leaves_names)
                        if len(tmp) == len(mrca_leaves_names):
                            mrcas_ids_to_clean.append(target_mrca)

                monophyletic_clades_leaves[mrca_id] = mrca_leaves_names
                monophyletic_clades[mrca_id] = mrca


    # remove MRCAs already contained in higher MRCAs
    mrcas_ids_to_clean = list(set(mrcas_ids_to_clean))
    for mrca_id in mrcas_ids_to_clean:
        del monophyletic_clades[mrca_id]
        del monophyletic_clades_leaves[mrca_id]


    # change mrca_ids to lower numbers
    final_monophyletic_clades = dict()
    final_cont = start_cont
    for mrca_id, mrca in monophyletic_clades.items():
        final_monophyletic_clades[f"mrca_{final_cont}"] = mrca
        if faces:
            mrca.add_features(mrca_id=f"mrca_{final_cont}")
            face = RectFace(5,5, "blue","blue", label=f"mrca_{final_cont}")
            mrca.add_face(face, 0, "float")
        final_cont +=1

    return t, final_monophyletic_clades, final_cont

def root_farthest_bacteria_hits(t, monophyletic_clades):
    # get the most distant node to the mrca
    rooted_trees = dict()
    for mrca_id, mrca in monophyletic_clades.items():
        bacteria_leaves_dists = {bacteria_leaf:mrca.get_distance(bacteria_leaf) for bacteria_leaf in t.search_nodes(superkingdom="Bacteria")}
        # check if there are still Bacteria after resolving polytomies
        if bacteria_leaves_dists:
            # sort by distance
            bacteria_leaves_dists = {node: distance for node, distance in sorted(bacteria_leaves_dists.items(), key=lambda item: item[1], reverse=True)}
            farthest_bact_leaf = list(bacteria_leaves_dists.keys())[0]
            # root using the farthest Bacteria
            t.set_outgroup(farthest_bact_leaf)
            #print(bacteria_leaves_dists)
            rooted_trees[farthest_bact_leaf.name] = t

    return rooted_trees

def process_MRCA(mrca):
    '''
    '''
    polytomy_seqs = ",".join([leaf.name for leaf in mrca.iter_leaves() if leaf.superkingdom == "bact_crass_polytomy"])
    n_seqs_mrca = len(mrca.get_leaf_names())

    sister = mrca.get_sisters()
    if len(sister) > 1:
        print("wuh?")
    n_seqs_sister = len(sister[0].get_leaf_names())
    bacteria_seqs_sister = ",".join([leaf.name for leaf in sister[0].iter_leaves() if leaf.superkingdom == "Bacteria"])
    sister_support = mrca.up.name

    return polytomy_seqs, n_seqs_mrca, n_seqs_sister, sister_support, bacteria_seqs_sister

def process_tree(crassvirales_taxa_df, taxa_groups, run_name, tree_file):

    # get CL_id
    cl_id = os.path.basename(tree_file).split("_ncbi")[0]

    # read summary file of the CL to know wheter it needs to be processed for Bacteria
    summary_dir = "/home/danielc/projects/Bas_phages/5_nr_screening/2_hits_summary"
    summary_file = f"{summary_dir}/{cl_id}.summary"
    ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0, low_memory=False)
    ncbi_summary_df = ncbi_summary_df[ncbi_summary_df.included]
    # check if there are Bacteria hits
    if "Bacteria" in ncbi_summary_df.superkingdom.to_list():
        # replace nan by "NA"
        ncbi_summary_df = ncbi_summary_df.fillna("NA")

        # parse vs2 results for the cl
        cl_vs2_prophages = parse_cl_vs2_results(cl_id)
        #print(cl_id, cl_vs2_prophages)


        # read the CL tree
        t = read_and_annotate_tree(tree_file, ncbi_summary_df, crassvirales_taxa_df, cl_vs2_prophages)
        # resolve Crassvirales-Bacteria polytomies
        t = resolve_crass_bacteria_polytomies(t)
        # identify Crassvirales MRCAs in the tree
        start_cont = 1
        t, mrcas, cont = identify_crassvirales_mrcas(t, taxa_groups, start_cont)
        # for all the MRCAs identified, get the farthest Bacteria sequence
        # and root with it. Return as many trees as different outgroups
        rooted_trees = root_farthest_bacteria_hits(t, mrcas)
        if rooted_trees:
            # for each of the rooted trees, indentify MRCAs again
            cont = 1
            header = ["mrca_id", "cl_id", "outgroup", "polytomies", "n_seqs_mrca", "n_seqs_sister", "sister_support", "bacteria_seqs_sister"]
            to_write = list()
            for out_seq, rooted_tree in rooted_trees.items():
                t, mrcas, cont = identify_crassvirales_mrcas(rooted_tree, taxa_groups, cont, faces=True)
                #print(mrcas)
                for mrca_id, mrca in mrcas.items():
                    polytomy_seqs, n_seqs_mrca, n_seqs_sister, sister_support, bacteria_seqs_sister = process_MRCA(mrca)

                    to_write.append([mrca_id,
                                    cl_id,
                                    out_seq,
                                    polytomy_seqs,
                                    n_seqs_mrca,
                                    n_seqs_sister,
                                    sister_support,
                                    bacteria_seqs_sister])

            # convert to df and write
            df = pd.DataFrame(to_write, columns=header)
            outfile = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{run_name}/0_mrcas/{cl_id}.txt"
            df.to_csv(outfile, header=True, index=False, sep="\t")


def main():

    args = parse_args()


    # create output dirs
    os.makedirs(f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{args.run_name}/0_mrcas", exist_ok=True)

    # get crassvirales taxonomy
    crassvirales_taxa_file = "/home/danielc/projects/Bas_phages/crass_protein_families_scripts/tree_functions/crassphages_taxonomy_terL_and_new.txt"
    crassvirales_taxa_df = get_crassvirales_taxonomy_df(crassvirales_taxa_file)

    # target taxa clades to look for in the trees
    taxa_groups = create_target_taxa_groups()
    taxa_groups = [["crassvirales"]] + taxa_groups # add crassvirales at the beginning

    # list trees
    trees_dir = "/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees"
    tree_files = glob.glob(f"{trees_dir}/*.nw")

    part = partial(process_tree, crassvirales_taxa_df, taxa_groups, args.run_name)
    pool = multiprocessing.Pool(processes=18)
    pool.map(part, tree_files)
    pool.close()
    pool.join()













if __name__ == "__main__":
    main()
