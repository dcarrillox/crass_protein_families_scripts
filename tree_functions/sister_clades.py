from ete3 import Tree, TreeStyle, TextFace, PhyloTree, NodeStyle, faces, AttrFace, CircleFace, RectFace
import pandas as pd
import random



def terl_genomes_seqs(terl_t):
    '''
    '''
    
    genomes_terl_seqs = {leaf.name.split("|")[0]:leaf.name for leaf in terl_t.iter_leaves()}
    
    return genomes_terl_seqs
        

def parse_sister_clades_table(sisters_clade_table):
    '''
    '''
    
    df = pd.read_csv(sisters_clade_table, header=0, sep="\t")
    
    return df


def merge_same_ancestor_cl(current_taxa, add_taxa):
    '''
    '''
    
    updated_taxa = dict()
    for add_taxa_name, add_taxa_value in add_taxa.items():
        if add_taxa_name in current_taxa:
            updated_taxa[add_taxa_name] = add_taxa_value + current_taxa[add_taxa_name]
        else:
            updated_taxa[add_taxa_name] = add_taxa_value
            
    for current_taxa_name, current_taxa_value in current_taxa.items():
        if current_taxa_name not in updated_taxa:
            updated_taxa[current_taxa_name] = current_taxa_value
            
    # get the average in the updated dict
    for taxa, value in updated_taxa.items():
        updated_taxa[taxa] = value/2
        
    return updated_taxa
        
    
    
    
    
def process_crass_clade(cl_id, clade_seqs_name, terl_genomes_seqs, terl_tree, taxa_sister):
    '''
    '''
    # check which sequences in the clade are Crassvirales and get their terl identifier
    clade_genomes = [name.split("|")[0] for name in clade_seqs_name]
    crass_terl_names = [terl_genomes_seqs[genome] for genome in clade_genomes if genome in terl_genomes_seqs]

    # identify the ancestor in the terl tree
    if len(crass_terl_names) == 1:
        clade_ancestor = terl_tree.search_nodes(name=crass_terl_names[0])[0]
    else:
        clade_ancestor = terl_tree.get_common_ancestor(crass_terl_names)
        
    
    # add taxa annotation
    taxas_split = taxa_sister.split(";")
    if len(taxas_split) > 1:
        taxas_values = {taxa.split(":")[0]:float(taxa.split(":")[1]) for taxa in taxas_split}
    else:
        taxas_values = {taxas_split[0].split(":")[0]:float(taxas_split[0].split(":")[1])}
        
        
    #if cl_id in clade_ancestor.phyla:
    #    print(f"{cl_id} already present in {clade_ancestor.name}")
    if cl_id in clade_ancestor.phyla:
        print(f"{cl_id} already present in {clade_ancestor.name}")
        #print(clade_ancestor.phyla)
        print(f"\t-Taxa before: {clade_ancestor.phyla[cl_id]}")
        print(f"\t-Taxa to add: {taxas_values}")
        updated_taxa = merge_same_ancestor_cl(clade_ancestor.phyla[cl_id], taxas_values)
        print(f"\t-Updated: {updated_taxa}")
        
        clade_ancestor.phyla[cl_id] = updated_taxa
        
        
            
              
    else:
        clade_ancestor.phyla[cl_id] = taxas_values
    #lade_ancestor.order[cl_id] = taxas_values
    
    
    
    #print(clade_ancestor.phyla)
    
    return terl_tree
    
    
def phyla_terl_nodes(terl_tree):
    '''
    '''
    
    for node in terl_tree.traverse():
        
        if node.phyla:

            taxa_df = pd.DataFrame.from_dict(node.phyla)
            taxa_df = taxa_df.fillna(0)
            taxa_df['mean'] = taxa_df.mean(axis=1)
            taxa_df = taxa_df.sort_values("mean", ascending=False)

            taxas_sorted = ";".join([f"{taxa}:{round(taxa_df.loc[taxa, 'mean'],2)}" for taxa in taxa_df.index])

            # remove mean column before reporting
            taxa_df.drop("mean", axis=1, inplace=True)
                
            node.add_features(ncl=taxa_df.shape[1],
                              taxa=taxas_sorted, 
                              cl=",".join(taxa_df.columns.tolist()))
            
    return terl_tree


def order_terl_nodes(terl_tree):
    '''
    '''
    
    for node in terl_tree.traverse():
        if node.order:

            taxa_df = pd.DataFrame.from_dict(node.order)
            taxa_df = taxa_df.fillna(0)
            taxa_df['mean'] = taxa_df.mean(axis=1)
            taxa_df = taxa_df.sort_values("mean", ascending=False)

            taxas_sorted = ";".join([f"{taxa}:{round(taxa_df.loc[taxa, 'mean'],2)}" for taxa in taxa_df.index])

            # remove mean column before reporting
            taxa_df.drop("mean", axis=1, inplace=True)
                
            node.add_features(ncl=taxa_df.shape[1],
                              taxa=taxas_sorted, 
                              cl=",".join(taxa_df.columns.tolist()))
            
    return terl_tree
            
            

def get_itol_annotation_phyla(terl_tree):
    '''
    '''
    
    itol_annot = list() 
    
    ## first I need to get all the possible phyla. 
    #all_phyla = list()
    #
    #for node in terl_tree.traverse():
    #    if not node.is_leaf():
    #        if node.phyla:
#
    #            taxa_df = pd.DataFrame.from_dict(node.phyla)
    #            all_phyla += taxa_df.index.tolist()
    #            
    #all_phyla = list(set(all_phyla))
    
    # Let's start easy, define the most important phylum by hand
    all_phyla = ["Actinobacteria", "Firmicutes", "Bacteroidetes", "Proteobacteria", "Spirochaetes"]
        
    for node in terl_tree.traverse():
        to_add = [node.name, "0" ]
        if node.phyla:

            taxa_df = pd.DataFrame.from_dict(node.phyla)
            # add number of cls as the pie diameter
            to_add.append(str(len(taxa_df.columns)))

            taxa_df = taxa_df.fillna(0)
            taxa_df['mean'] = taxa_df.mean(axis=1)
            taxa_df = taxa_df.sort_values("mean", ascending=False)
                
                
            for phylum in all_phyla:
                if phylum in taxa_df.index:
                    to_add.append(str(round(taxa_df.loc[phylum, "mean"],2)))
                else:
                    to_add.append("0")
                        
                        
            others = 0
            for phylum in taxa_df.index:
                if phylum not in all_phyla:
                    others += taxa_df.loc[phylum, "mean"]
            to_add.append(str(round(others, 2)))

            itol_annot.append(to_add)
            
            
    # remove annot for nodes without CLs
    itol_annot = [node_annot for node_annot in itol_annot if len(node_annot) > 2]
    all_phyla.append("Others")
    
    return itol_annot, all_phyla
    

def get_itol_annotation_order(terl_tree):
    '''
    '''
    
    itol_annot = list() 
    
    # start by getting all the possible taxa names
    
    ## first I need to get all the possible phyla. 
    all_orders = list()
    
    for node in terl_tree.traverse():
        if node.order:
            taxa_df = pd.DataFrame.from_dict(node.order)
                
            taxa_df = taxa_df.fillna(0)
            taxa_df['mean'] = taxa_df.mean(axis=1)
            taxa_df = taxa_df.sort_values("mean", ascending=False)
                
            taxa_df = taxa_df[taxa_df["mean"] > 0.05]
            all_orders += taxa_df.index.tolist()
                
    all_orders = list(set(all_orders))
    
    # Iterate the nodes again, this time assigning values to the orders that were higher 
    # than 5% in any CL
   
        
    for node in terl_tree.traverse():
        to_add = [node.name, "0" ]
        if node.order:

            taxa_df = pd.DataFrame.from_dict(node.order)
            # add number of cls as the pie diameter
            to_add.append(str(len(taxa_df.columns)))

            taxa_df = taxa_df.fillna(0)
            taxa_df['mean'] = taxa_df.mean(axis=1)
            taxa_df = taxa_df.sort_values("mean", ascending=False)
                
                
            for order in all_orders:
                if order in taxa_df.index:
                    to_add.append(str(round(taxa_df.loc[order, "mean"],2)))
                else:
                    to_add.append("0")
                        
                        
            others = 0
            for order in taxa_df.index:
                if order not in all_orders:
                    others += taxa_df.loc[order, "mean"]
            to_add.append(str(round(others, 2)))

            itol_annot.append(to_add)
            
            
    # remove annot for nodes without CLs
    itol_annot = [node_annot for node_annot in itol_annot if len(node_annot) > 2]
    all_orders.append("Others")
    
    return itol_annot, all_orders
    

def generate_random_hex_color():
    '''
    '''
    hex_color = "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])
    
    return hex_color
    

def write_itol_taxa_annot_terl(taxa_annot_df, terl_tree, outfile):
    '''
    '''
    
    families_colors = {"Intestiviridae":"#EE3B3B",
                      "Crevaviridae":"#EE9A00",
                      "Suoliviridae":"#4169E1", 
                      "Steigviridae":"#00CED1",
                      "Tinaiviridae":"#CD2990",
                      "Jelitoviridae":"#006400"
                     }
       
    
    to_write = ["TREE_COLORS\nSEPARATOR TAB\nDATA"]
    #leaves_names = [leaf.name for leaf in terl_tree.get_leaves()]
    #for leaf_name in leaves_names:
    #    genome = leaf_name.split("|")[0]
    #    family = taxa_annot_df.loc[genome, "family"]
    #    if family in colors:
    #        to_write.append(f"{leaf_name}\trange\t{colors[family]}\t{family}")
    #
   
    
    # identify the clades of each family in the tree
    for family, color in families_colors.items():
        print(family)
        fam_leaves = terl_tree.search_nodes(family=family)
        fam_lca = terl_tree.get_common_ancestor(fam_leaves)
        to_write.append(f"{fam_lca.name}\tclade\t{families_colors[family]}\tnormal\t2")
        
    with open(outfile, "w") as fout:
        for line in to_write:
            fout.write(line + "\n")

       
    
    
def write_itol_dataset(itol_annot, labels, template, outfile):
    '''
    '''
  
    # read template
    lines = [line.strip() for line in open(template).readlines()]
    
    # replace name, labels
    for i, line in enumerate(lines):
        if line == "SEPARATOR COMMA":
            lines[i] = "#SEPARATOR COMMA"
            
        if line == "#SEPARATOR TAB":
            lines[i] = "SEPARATOR\tTAB"
        
        if line == "FIELD_LABELS,f1,f2,f3":
            lines[i] = "FIELD_LABELS\t" + "\t".join(labels) 
        
        if line == "FIELD_COLORS,#ff0000,#00ff00,#0000ff":
            lines[i] = "FIELD_COLORS\t" + "\t".join([generate_random_hex_color() for i in range(len(labels))])
                
        if line == "#LEGEND_TITLE,Dataset legend":
            lines[i] = "LEGEND_TITLE\tPhyla"
            
        if line == "COLOR,#ff0000":
            lines[i] = "COLOR\t#ff0000"
            
        if line == "DATASET_LABEL,example multi bar chart":
            lines[i] == "DATASET_LABEL\tPhyla"
    
    
    # add annot lines
    lines += itol_annot
    
    
    # write to file
    with open(outfile, "w") as fout:
        for line in lines:
            #rint(line, type(line))
          
            if isinstance(line, str):
                fout.write(line + "\n")
            else:
                line = "\t".join(line)
                fout.write(line + "\n")
            