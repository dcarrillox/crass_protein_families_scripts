import os, glob

# list hmm files
hmm_files = glob.glob("*.hmm")

outdirs = ["hmmsearch1", "hmmsearch2", "hmmsearch3", "hmmsearch4"] * 5000

i = 0
for hmm_file in hmm_files:
    os.system(f"cp {hmm_file} {outdirs[i]}")
    i += 1
