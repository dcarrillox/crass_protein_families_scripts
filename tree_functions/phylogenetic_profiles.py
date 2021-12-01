from ete3 import Tree, TreeStyle, TextFace, PhyloTree, NodeStyle, faces, AttrFace, CircleFace, RectFace
import pandas as pd



def cls_genomes(cls_proteins_file):
    '''
    It reads the clusters_table generated with 1_mmseqs2tsv_to_faa.py and saves in
    a dict k=cl_id  v=[genome1, genome2, genome3]
    '''

    lines = [line.strip().split("\t") for line in open(cls_proteins_file).readlines()]

    cls_genomes = {line[0]:list() for line in lines}
    for line in lines:
        cls_genomes[line[0]] += [line[1].split("|")[0]]

    # remove redundancy
    for cl, genomes in cls_genomes.items():
        cls_genomes[cl] = list(set(genomes))

    return cls_genomes


def color_profiles_red_palette(value):
    '''
    '''
    if value <= 1:
        color = "#a70000"
    if value < 0.8:
        color = "#ff0000"
    if value < 0.6:
        color = "#ff5252"
    if value < 0.4:
        color = "#ff7b7b"
    if value < 0.2:
        color = "#ffbaba"
    if value < 0.1:
        color = "#FFD4D4"
    if value == 0:
        color = "#FFFFFF"

    return color


def add_profiles_representatives(tree, genomes_reprs, votus, target_cls, cls_genomes):
    '''
    In a terl tree made with the representatives of the vOTUs, it adds RectFace faces
    to the tree colored according to the presence of the target_cls in the whole vOTU.
    '''

    # red color palette from https://www.color-hex.com/color-palette/5634

    column_face = 1
    for cl in target_cls:
        # for the genomes where the cl is present, get their representative genome in the vOTU
        reprs = [genomes_reprs[genome] for genome in cls_genomes[cl]]
        reprs_uniq = list(set(reprs))

        # iterate the reprs while assigning a color to them based on the n_genomes with the cl
        reprs_color = dict()
        reprs_metadata = dict()
        for repr in reprs_uniq:
            n_genomes = reprs.count(repr)
            fraction = float(n_genomes/votus[repr]["size"])
            reprs_color[repr] = {"color":color_profiles_red_palette(fraction), "n_genomes":str(n_genomes)}


        # iterate the tree while adding the profiles
        for leaf in tree.iter_leaves():
            if leaf.genome in reprs_color:
                color = reprs_color[leaf.genome]["color"]
                label = {"text":reprs_color[leaf.genome]["n_genomes"], "color":"black", "fontsize": 7}
                face = RectFace(10, 10, "white", color, label=label)
            else:
                color = "#C0C0C0"
                face = RectFace(10, 10, "white", color)


            leaf.add_features(nvotu=votus[leaf.genome]["size"])
            leaf.add_face(face, column_face, "aligned")

        column_face += 1

    return tree
