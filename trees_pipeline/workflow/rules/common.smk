import glob, os

configfile: "config/config.yaml"

# get the clusters_id of the .faa files in the input folder
input_faa_files = glob.glob(f"{config['input_dir']}/*.faa")
cls_files = {os.path.basename(faa_file).replace(".faa", ""):faa_file for faa_file in input_faa_files}

CLUSTERS = list(cls_files.keys())


def get_input_faa(wildcards):
    return cls_files[wildcards.cl]
