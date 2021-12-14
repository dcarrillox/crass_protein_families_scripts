import os, glob
import pandas as pd
from itertools import combinations
from collections import Counter
from ete3 import Tree, TreeStyle, TextFace, PhyloTree, NodeStyle, faces, AttrFace, CircleFace, RectFace
from Bio import SeqIO

def parse_ncbi_summary_table(cl_id, summaries_dir):
    '''
    For a given CL id, parses the summary table with the NCBI hits for the CL
    and returns a DataFrame with k=ncbi_prot_id
    v={source=ncbi, crass=True/False, superkingdom: , phylum: , class: , ...}
    '''
    summary_file = f"{summaries_dir}/{cl_id}.summary"

    ncbi_summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
    # remove included=False rows
    ncbi_summary_df = ncbi_summary_df[ncbi_summary_df.included]
    # replace nan by "NA"
    ncbi_summary_df = ncbi_summary_df.fillna("NA")
    return ncbi_summary_df

def read_and_annotate_tree(cl_id, trees_dir, ncbi_summary_df, crassvirales_taxa_df):
    '''
    Reads the tree of a CL with NCBI and Crassvirales sequences and assigns
    source and taxonomy to them
    '''

    # read tree
    tree_file = f"{trees_dir}/{cl_id}_ncbi_trimmed.nw"
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

def add_phylum_percent_inner_nodes(t):
    '''

    '''
    for node in t.traverse():
        # process only inner nodes
        if not node.is_leaf():
            # get all phylums in the childs
            phylums = [child.phylum for child in node.iter_leaves()]
            # get the percentage of each phylum
            counter = Counter(phylums)
            total = sum(counter.values())
            # store each phylum along with its percentage and count
            phylums_percent = {phylum:[round(count/total, 3), count] for phylum, count in counter.items()}

            node.add_features(phylum_perc=phylums_percent)

    return t

def collapse_nodes_phylum_perc(t, perc_cutoff):
    '''

    '''
    for node in t.traverse():
        if not node.is_leaf():
            # check if a phylum accounts for more than 90% of the leaves
            major_phylum = [(phylum, values) for phylum, values in node.phylum_perc.items() if values[0] >= perc_cutoff]
            # if so, detach childs
            if major_phylum:
                # set the phylum for the node
                node.add_features(phylum=major_phylum[0][0])
                phylum = major_phylum[0][0]
                perc  = major_phylum[0][1][0]
                count = major_phylum[0][1][1]
                node.name = f"{phylum}, {perc}, {count}"

                # get the name of the childs, add it as a feature
                child_leaves = [child.name for child in node.iter_descendants()]
                node.add_features(childs=",".join(child_leaves))
                # collapse nodes below this node
                for child in node.iter_descendants():
                    child_leaves.append
                    child.detach()
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

def get_interesting_cls(trees_dir, summaries_dir, cls_faa_dir, crassvirales_tax_df):
    '''

    '''
    # get the CLs for which I have a non-empty tree
    ncbi_cls = [os.path.basename(file).split("_ncbi")[0] for file in glob.glob(f"{trees_dir}/*.nw") if os.path.getsize(file) != 0]

    interesting_cls = list()
    # for these CLs, parse the NCBI summary and Crassvirales taxa to know which
    # ones contain more than one family
    for cl in ncbi_cls:
        summary_file = f"{summaries_dir}/{cl}.summary"
        summary_df = pd.read_csv(summary_file, sep="\t", header=0, index_col=0)
        summary_df = summary_df[summary_df.included]
        #print(cl, summary_df)

        if "Bacteria" in summary_df.superkingdom.values:
            # check if there are several families of Crassvirales in the cl
            cl_faa_file = f"{cls_faa_dir}/{cl}.faa"
            records = SeqIO.parse(cl_faa_file, "fasta")
            genomes = [record.name.split("|")[0] for record in records]
            crassv_families = list(set([crassvirales_tax_df.loc[genome, "family"] for genome in genomes]))
            if "new" in crassv_families:
                crassv_families.remove("new")

            if len(crassv_families) > 1:
                print(cl, "Bacteria", summary_df[summary_df.superkingdom == 'Bacteria'].shape[0])




# def get_taxa_names_rank(t, rank):
#     '''
#
#     '''
#
#     if rank == "phylum":
#         taxa_names_rank = list(set([leaf.phylum for leaf in t.iter_leaves()]))
#
#     return taxa_names_rank


# def collapse_nodes_XXpercent(t, rank, percent):
#     '''
#     Collapses nodes at the specified rank. First it gets all the taxas in the tree
#     for the rank specified. Then creates all the possible combinations of names
#     and identifies the monophyletic ranks with these names. If one of the names makes
#     up to a certain percentage of the clade, the node is collapsed with that name.
#     '''
#
#     # get all the names for the rank in the tree
#     taxa_names_rank = get_taxa_names_rank(t, rank)
#
#     taxa_names_combinations = list(combinations(taxa_names_rank, 2))
