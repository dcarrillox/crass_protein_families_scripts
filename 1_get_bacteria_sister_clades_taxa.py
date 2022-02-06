'''

'''

import os, glob, argparse, multiprocessing, time, itertools
from functools import partial
import pandas as pd
import numpy as np
from collections import Counter
from ete3 import Tree, RectFace


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-td', '--trees_dir',
                               dest='trees_dir',
                               required=True,
                               help='directory with the .nw trees'
                               )
    requiredArgs.add_argument('-s', '--summary_dir',
                               dest='summary_dir',
                               required=True,
                               help='directory to save the .summary files'
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

def color_format_tree(t):
    '''

    '''

    crassvirales_colors = {'Crevaviridae': 'red',
                           'Intestiviridae': 'cyan',
                           'Jelitoviridae': 'green',
                           'Steigviridae': 'orange',
                           'Suoliviridae': 'violet',
                           'Tinaiviridae': 'brown'}

    for node in t.traverse():
        # check if it is a Crassvirales genome
        if node.is_leaf():
            if node.phylum in crassvirales_colors:
                node.img_style["bgcolor"] = crassvirales_colors[node.phylum]
            else:
                node.img_style["size"] = 0

            # color Bacteria as yellow
            if node.superkingdom == "Bacteria":
                node.img_style["bgcolor"] = "yellow"

    return t

def read_and_annotate_tree(tree_file, ncbi_summary_df, crassvirales_taxa_df):
    '''
    Reads the tree of a CL with NCBI and Crassvirales sequences and assigns
    source and taxonomy to them
    '''

    colors = {}

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

def identify_crassvirales_mrcas(t, cl_id, taxa_groups, start_cont, faces=False):
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
            face = RectFace(5,5, "red","red", label=f"mrca_{final_cont}")
            mrca.add_face(face, 0, "float")
        final_cont +=1

    if faces:
        pdf_outfile = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/4_MRCAs_taxa/{cl_id}.pdf"
        t.render(pdf_outfile)

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

def get_taxa_sister(sister_clade):
    '''
    '''
    superkingdoms = [leaf.superkingdom for leaf in sister_clade.iter_leaves()]

    # get the percentage of each phylum
    counter = Counter(superkingdoms)
    total = sum(counter.values())
    # store each phylum along with its percentage and count
    superkingdoms_percent = [(superkingdom, round(count/total, 3)) for superkingdom, count in counter.items()]
    superkingdoms_percent = sorted(superkingdoms_percent, key=lambda superkingdom: superkingdom[1], reverse=True )
    superkingdoms_percent = ";".join([f"{superkingdom[0]}:{superkingdom[1]}" for superkingdom in superkingdoms_percent])


    # store the phyla for the Bacteria hits
    phyla = [leaf.phylum for leaf in sister_clade.iter_leaves() if leaf.superkingdom == "Bacteria"]
    counter = Counter(phyla)
    total = sum(counter.values())
    # store each phylum along with its percentage and count
    phyla_percent = [(phylum, round(count/total, 3)) for phylum, count in counter.items()]
    phyla_percent = sorted(phyla_percent, key=lambda phylum: phylum[1], reverse=True )
    phyla_percent = ";".join([f"{phylum[0]}:{phylum[1]}" for phylum in phyla_percent])


    # store order for the Bacteria hits
    order = [f"{leaf.phylum}_{leaf.order}" for leaf in sister_clade.iter_leaves() if leaf.superkingdom == "Bacteria"]
    counter = Counter(order)
    total = sum(counter.values())
    # store each phylum along with its percentage and count
    order_percent = [(order, round(count/total, 3)) for order, count in counter.items()]
    order_percent = sorted(order_percent, key=lambda order: order[1], reverse=True )
    order_percent = ";".join([f"{order[0]}:{order[1]}" for order in order_percent])

    return superkingdoms_percent, phyla_percent, order_percent

def process_MRCA(mrca):
    '''
    '''
    polytomy_seqs = ",".join([leaf.name for leaf in mrca.iter_leaves() if leaf.superkingdom == "bact_crass_polytomy"])
    seqs_mrca = ",".join(mrca.get_leaf_names())
    n_seqs_mrca = len(mrca.get_leaf_names())

    sister = mrca.get_sisters()
    if len(sister) > 1:
        print("wuh?")
    seqs_sister = ",".join(sister[0].get_leaf_names())
    n_seqs_sister = len(sister[0].get_leaf_names())
    n_seqs_bacteria = len([leaf.name for leaf in sister[0].iter_leaves() if leaf.superkingdom == "Bacteria"])
    bacteria_seqs_sister = ",".join([leaf.name for leaf in sister[0].iter_leaves() if leaf.superkingdom == "Bacteria"])
    sister_support = mrca.up.name

    superkingdoms_sister, bact_phyla_sister, bact_order_sister = get_taxa_sister(sister[0])

    return polytomy_seqs, n_seqs_mrca, seqs_mrca, sister_support, n_seqs_sister, seqs_sister, n_seqs_bacteria, bacteria_seqs_sister, superkingdoms_sister, bact_phyla_sister, bact_order_sister

def process_tree(summary_dir, crassvirales_taxa_df, taxa_groups, tree_file):

    # get CL_id
    cl_id = os.path.basename(tree_file).split("_ncbi")[0]

    # read summary file of the CL to know wheter it needs to be processed for Bacteria
    summary_file = f"{summary_dir}/{cl_id}.summary"
    ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0, low_memory=False)
    ncbi_summary_df = ncbi_summary_df[ncbi_summary_df.included]
    # check if there are Bacteria hits
    if "Bacteria" in ncbi_summary_df.superkingdom.to_list():
        # replace nan by "NA"
        ncbi_summary_df = ncbi_summary_df.fillna("NA")
        #print(ncbi_summary_df)

        # read the CL tree
        t = read_and_annotate_tree(tree_file, ncbi_summary_df, crassvirales_taxa_df)
        # add colors
        t = color_format_tree(t)
        # resolve Crassvirales-Bacteria polytomies
        t = resolve_crass_bacteria_polytomies(t)
        # identify Crassvirales MRCAs in the tree
        start_cont = 1
        t, mrcas, cont = identify_crassvirales_mrcas(t, cl_id, taxa_groups, start_cont)
        # for all the MRCAs identified, get the farthest Bacteria sequence
        # and root with it. Return as many trees as different outgroups
        rooted_trees = root_farthest_bacteria_hits(t, mrcas)
        if rooted_trees:
            # for each of the rooted trees, indentify MRCAs again
            cont = 1
            header = ["mrca_id", "cl_id", "outgroup", "polytomies", "n_seqs_mrca", "seqs_mrca", "sister_support",
                      "n_seqs_sister", "seqs_sister", "n_seqs_bacteria", "bacteria_seqs_sister",
                      "superkingdoms_sister", "bact_phyla_sister", "bact_order_sister"]
            to_write = list()
            for out_seq, rooted_tree in rooted_trees.items():
                t, mrcas, cont = identify_crassvirales_mrcas(rooted_tree, cl_id, taxa_groups, cont, faces=True)
                #print(mrcas)
                for mrca_id, mrca in mrcas.items():
                    polytomy_seqs, n_seqs_mrca, seqs_mrca, sister_support, n_seqs_sister, seqs_sister, n_seqs_bacteria, bacteria_seqs_sister, superkingdoms_sister, bact_phyla_sister, bact_order_sister = process_MRCA(mrca)

                    to_write.append([mrca_id,
                                    cl_id,
                                    out_seq,
                                    polytomy_seqs,
                                    n_seqs_mrca,
                                    seqs_mrca,
                                    sister_support,
                                    n_seqs_sister,
                                    seqs_sister,
                                    n_seqs_bacteria,
                                    bacteria_seqs_sister,
                                    superkingdoms_sister,
                                    bact_phyla_sister,
                                    bact_order_sister])

            # convert to df and write
            df = pd.DataFrame(to_write, columns=header)
            outfile = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/4_MRCAs_taxa/{cl_id}.txt"
            df.to_csv(outfile, header=True, index=False, sep="\t")


def main():

    args = parse_args()


    # create output dirs
    #os.makedirs(f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/3_bacterial/{args.run_name}/0_mrcas", exist_ok=True)

    # get crassvirales taxonomy
    crassvirales_taxa_file = "/home/danielc/projects/Bas_phages/crass_protein_families_scripts/tree_functions/crassphages_taxonomy_terL_and_new.txt"
    crassvirales_taxa_df = get_crassvirales_taxonomy_df(crassvirales_taxa_file)

    # target taxa clades to look for in the trees
    taxa_groups = create_target_taxa_groups()
    taxa_groups = [["crassvirales"]] + taxa_groups # add crassvirales at the beginning

    # list trees
    tree_files = glob.glob(f"{args.trees_dir}/*.nw")

    part = partial(process_tree, args.summary_dir, crassvirales_taxa_df, taxa_groups)
    pool = multiprocessing.Pool(processes=18)
    pool.map(part, tree_files)
    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
