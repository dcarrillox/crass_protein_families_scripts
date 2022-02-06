import os, glob
import pandas as pd
from itertools import combinations
from collections import Counter
from ete3 import Tree, TreeStyle, TextFace, PhyloTree, NodeStyle, faces, AttrFace, CircleFace, RectFace
from Bio import SeqIO


def parse_vs2_annots(dir):
    '''
    '''

    # list all the classification tables
    classification_files = glob.glob(f"{dir}/round*/3_vs2/2_classified/*.prophages")

    # parse them with pandas
    all_dfs = list()
    for file in classification_files:
        df = pd.read_csv(file, sep="\t", header=0, low_memory=False)
        # keep only analyzed proteins
        df = df[(df["ipg_id"] != 0) & (df.ipg_id.notnull()) & (df["prophage"] != "0.0")]
        
        if not df.empty:
            all_dfs.append(df)

    proteins_df = pd.concat(all_dfs)
    proteins_df.drop_duplicates(ignore_index=True, inplace=True)

    proteins_df = proteins_df.rename(columns={"Unnamed: 0": "prot_id"})
    proteins_df.set_index("prot_id", inplace=True)
    proteins_df = proteins_df[~proteins_df.index.duplicated(keep='first')]

    
    # parse contigs information (full/partial)
    classification_files = glob.glob(f"{dir}/round*/3_vs2/1_vs2_results/final-viral-boundary*.tsv")
    types_list = list()
    for file in classification_files:
        df = pd.read_csv(file, sep="\t", header=0)
        # split seqname_new column by ||
        contigs_types = [[seq.split("||")[0], seq.split("||")[1]] for seq in df["seqname_new"].tolist()]
        for contig_type in contigs_types:
            if "partial" in contig_type[1]:
                contig_type[1] = "partial"

        types_list += contigs_types

    types_df = pd.DataFrame(types_list, columns=["contig", "type"])
    types_df.drop_duplicates(ignore_index=True, inplace=True)
    types_df.set_index("contig", inplace=True)
    
    # add contig type to protein annotation
    def add_contig_type(protein_line):

        if protein_line["prophage"]:
            #print(protein_line)
            contigs = protein_line["prophage_contigs"].split(",")
            #print(contigs)
            types = list({types_df.loc[contig, "type"] for contig in contigs if contig in types_df.index})

            if not types:
                return "full"

            else:
                if len(types) > 1:
                    return "prophage_partial"
                else:
                    return f"prophage_{types[0]}"

        else:
            if protein_line["n_short"] == protein_line["n_contigs"]:
                return "short"
            else:
                return "bacteria"
                



    proteins_df["contig_type"] = proteins_df.apply(lambda x: add_contig_type(x), axis=1)


    return proteins_df

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
    print(tree_file)
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


def annotate_inner_nodes_superkingdom(t):
    '''
    '''
    for node in t.traverse():
        if not node.is_leaf():
            # get all phylums in the childs
            superkingdoms = [child.superkingdom for child in node.iter_leaves()]

            # get the percentage of each phylum
            counter = Counter(superkingdoms)
            total = sum(counter.values())
            # store each phylum along with its percentage and count
            superkingdoms_percent = {superkingdom:[round(count/total, 3), count] for superkingdom, count in counter.items()}

            node.add_features(kingdom_perc=superkingdoms_percent)
            
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

    crassvirales_colors = {"Intestiviridae":"#EE3B3B",
                      "Crevaviridae":"#EE9A00",
                      "Suoliviridae":"#4169E1", 
                      "Steigviridae":"#00CED1",
                      "Tinaiviridae":"#CD2990",
                      "Jelitoviridae":"#006400"
                     }

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
