import glob
import os
import multiprocessing
import time


def timeit(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end


# count the number of sequences in each CL_xx.faa and sort them from max to min
os.system("grep -c '>' *.faa > n-seqs.txt")
lines = [line.strip().split(".faa:") for line in open("n-seqs.txt").readlines()]
order = sorted(lines, key=lambda CL: int(CL[1]), reverse=True)

# define out dir for the alignments
out_dir = "/home/danielc/projects/Bas_phages/3_Broccoli/4_easy_CLs/1_easy_CLs_mafft"
threads = 90

files = glob.glob("*faa")

for CL in order:
    start = time.time()
    print(CL)

    ## ALWAYS RUN MAFFT-EINSI
    os.system(f"mafft-einsi --thread {threads} {CL[0]}.faa > {out_dir}/{CL[0]}.mafft 2> {out_dir}/{CL[0]}.log")


    # ## RUN MAFFT-FFTSNI 1000 IF LOT OF SEQUENCES
    # if int(CL[1]) < 1000:
    #     os.system(f"mafft-einsi --thread {threads} {CL[0]}.faa > {out_dir}/{CL[0]}.mafft 2> {out_dir}/{CL[0]}.log")
    # else:
    #     print("FFT-NS-I 200iter")
    #     os.system(f"mafft-fftnsi --maxiterate 200 --thread {threads} {CL[0]}.faa > {out_dir}/{CL[0]}.mafft 2> {out_dir}/{CL[0]}.log")



    timeit(start)
    print()
