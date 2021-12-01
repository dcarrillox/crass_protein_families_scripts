'''
Separates the CLs.faa in different folders to run mafft in different mutants
simultaneously.

There are some CLs that I don't want to align yet. They are in the list "not_align"
'''

import glob
import os

not_align = []

folders = ["mafft1", "mafft2"] * 5000

os.system("grep -c '>' *faa > n-seqs.txt")
lines = [line.strip().split(".faa:") for line in open("n-seqs.txt").readlines()]

order = sorted(lines, key=lambda CL: int(CL[1]), reverse=True)

#print(order)

i = 0
for CL in order:
    if CL[0].split("_")[1] not in not_align:
        os.system(f"cp {CL[0]}.faa {folders[i]}")
        i += 1
    else:
        print(CL)
