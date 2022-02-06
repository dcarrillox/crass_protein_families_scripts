'''
For a given CL, parses the tree to get the Crassvirales monophyletic clades and
the Bacteria (if any) in the sister clade. For each of the Crassvirales clades,
root the tree with the most distant Bacteria sequence.
'''

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import argparse
import warnings
from ete3 import Tree, RectFace

import glob, os

def get_taxonomy_df(taxa_file):
    '''
    Parses the taxonomy and returns a dataframe, index=genomes, 3 columns (family, subfamily, genus)
    '''
    #"/home/danielc/software/github/old_june/tree_functions/files/crassphages_taxonomy_terL_and_new.txt"
    df = pd.read_csv(taxa_file, sep="\t", index_col=0, names=["family","subfamily","genus"])
    df['superkingdom'] = np.where(df['family'].isin(["new", "crass_env"]), "new", "crassvirales")
    return df

def read_and_annotate_tree(tree_file, ncbi_summary_df, crassvirales_taxa_df):
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

    return t

def identify_crassvirales_mrcas(t, start_cont, faces=False):
    '''

    '''
    taxa_groups = [["crassvirales","NA", "new", "Viruses", "bact_crass_polytomy"],
               ["crassvirales", "new", "NA", "Viruses"],
               ["crassvirales", "new", "NA", "bact_crass_polytomy"],
               ["crassvirales", "new", "Viruses", "bact_crass_polytomy"],
               ["crassvirales", "NA", "Viruses", "bact_crass_polytomy"],
               ["crassvirales", "new", "NA"],
               ["crassvirales", "new", "Viruses"],
               ["crassvirales", "new", "bact_crass_polytomy"],
               ["crassvirales", "NA", "Viruses"],
               ["crassvirales", "NA", "bact_crass_polytomy"],
               ["crassvirales", "Viruses", "bact_crass_polytomy"],
               ["crassvirales", "new"],
               ["crassvirales", "NA"],
               ["crassvirales", "Viruses"],
               ["crassvirales", "bact_crass_polytomy"],
               ["crassvirales"]]

    monophyletic_clades_leaves = dict()
    monophyletic_clades = dict()
    mrcas_ids_to_clean = list()
    cont = 0
    for taxa_group in reversed(taxa_groups):
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
    mrcas_ids_to_clean   = list(set(mrcas_ids_to_clean))
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




def main():

    crassvirales_taxa_file = snakemake.config["crassvirales_taxa_file"]
    crassvirales_taxa_df = get_taxonomy_df(crassvirales_taxa_file)

    cl_id = os.path.basename(snakemake.input[0]).split("_ncbi_trimmed.nw")[0]

    summary_file = f"{snakemake.config['summaries_dir']}/{cl_id}.summary"

    if os.path.getsize(snakemake.input[0]) != 0:
        # parse the NCBI summary
        ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
        ncbi_summary_df = ncbi_summary_df[ncbi_summary_df.included]
        # check if there are Bacteria hits
        if "Bacteria" in ncbi_summary_df.superkingdom.values:
            # replace nan by "NA"
            ncbi_summary_df = ncbi_summary_df.fillna("NA")

            # read the CL tree
            t = read_and_annotate_tree(snakemake.input[0],ncbi_summary_df, crassvirales_taxa_df)
            # resolve Crassvirales-Bacteria polytomies
            t = resolve_crass_bacteria_polytomies(t)
            # identify Crassvirales MRCAs in the tree
            start_cont = 1
            t, mrcas, cont = identify_crassvirales_mrcas(t, start_cont)
            # for all the MRCAs identified, get the farthest Bacteria sequence
            # and root with it. Return as many trees as different outgroups
            rooted_trees = root_farthest_bacteria_hits(t, mrcas)
            if rooted_trees:
    	        # for each of the rooted trees, indentify MRCAs again
    	        cont = 1
    	        header = ["mrca_id", "cl_id", "outgroup", "polytomies", "n_seqs_mrca", "n_seqs_sister", "sister_support", "bacteria_seqs_sister"]
    	        to_write = list()
    	        for out_seq, rooted_tree in rooted_trees.items():
    	            t, mrcas, cont = identify_crassvirales_mrcas(rooted_tree, cont, faces=True)
    	            print(mrcas)
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
    	        df.to_csv(snakemake.output[0], header=True, index=False, sep="\t")
            else:
                os.system(f"touch {snakemake.output[0]}")

        else:
            os.system(f"touch {snakemake.output[0]}")


    else:
        os.system(f"touch {snakemake.output[0]}")

if __name__ == "__main__":
    main()
