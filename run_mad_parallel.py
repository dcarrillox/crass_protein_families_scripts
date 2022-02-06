import os, glob
import multiprocessing


# list the trees
trees_files = glob.glob("/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees/*.nw")


def root_tree(tree_file):
	
	mad_tree = f"{tree_file}.rooted"

	out_file = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_mad/{os.path.basename(mad_tree)}"
	log = f"/home/danielc/projects/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_mad/{os.path.basename(mad_tree)}.log"
	
	os.system(f"python /home/danielc/projects/Bas_phages/mad/mad.py {tree_file} > {log}")
	os.system(f"mv {mad_tree} {out_file}")




pool = multiprocessing.Pool(processes=15)
pool.map(root_tree, trees_files)
pool.close()
pool.join()

