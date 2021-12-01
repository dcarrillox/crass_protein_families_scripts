from ete3 import Tree, TreeStyle, TextFace, PhyloTree, NodeStyle, faces, AttrFace, CircleFace, RectFace
import pandas as pd


def get_taxonomy_df(taxa_file):
    '''
    Parses the taxonomy and returns a dataframe, index=genomes, 3 columns (family, subfamily, genus)
    '''
    #"/home/danielc/software/github/old_june/tree_functions/files/crassphages_taxonomy_terL_and_new.txt"
    df = pd.read_csv(taxa_file, sep="\t", index_col=0, names=["family","subfamily","genus"])
    return df


def parse_terl_tree(tree_file, tax_df, names=False, circular=False):
    '''
    Parses the tree_file to a Tree object. Tree is rooted with the outgroup,
    branches are colored based on family. Branch support values are grabbed from
    the name of the inner node.
    '''
    t = Tree(tree_file, format=1)

    for node in t.traverse():
        if node.is_leaf():
            genome = node.name.split("|")[0]
            node.add_features(family=tax_df.loc[genome, "family"],
                              subfamily=tax_df.loc[genome, "subfamily"],
                              genus=tax_df.loc[genome, "genus"],
                              genome=genome)
        else:
            # branch support
            if node.name != "":
                node.support = round(float(node.name), 2)

        node.img_style["size"] = 0

    outgs_leaves = t.search_nodes(family="outgroup")
    outgs_lca = t.get_common_ancestor(outgs_leaves)
    t.set_outgroup(outgs_lca)

    families_colors = {'Crevaviridae': 'red',
                       'Intestiviridae': 'cyan',
                       'Jelitoviridae': 'green',
                       'Steigviridae': 'orange',
                       'Suoliviridae': 'violet',
                       'Tinaiviridae': 'brown'}

    for family, color in families_colors.items():
        fam_leaves = t.search_nodes(family=family)
        fam_lca = t.get_common_ancestor(fam_leaves)
        fam_lca.name = family
        for node in fam_lca.traverse():
            #node.detach()
            node.img_style['hz_line_color'] = families_colors[family]
            node.img_style['vt_line_color'] = families_colors[family]
            node.img_style['hz_line_width'] = 2
            node.img_style['vt_line_width'] = 2

    ts = TreeStyle()
    ts.show_branch_support = False
    if circular:
        ts.mode = "c"

    if not names:
        ts.show_leaf_name = False

        for leaf in t.iter_leaves():
            face = TextFace(f'{tax_df.loc[leaf.genome, "subfamily"]},{tax_df.loc[leaf.genome, "genus"]}')
            leaf.add_face(face, column=0, position = "branch-right")


    return t, ts
