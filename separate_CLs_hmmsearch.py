import os, glob

'''
Changed to run with OGs from v7
'''

# list hmm files
hmm_files = glob.glob("*.hmm")

outdirs = ["hmmsearch1", "hmmsearch2", "hmmsearch3", "hmmsearch4", "hmmsearch5"] * 5000

path = "/home/danielc/projects/Bas_phages/4_scan_databases/0_hmm_profiles"

i = 0
for hmm_file in hmm_files:
    os.system(f"cp {hmm_file} {path}/{outdirs[i]}")
    i += 1
