'''

'''

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
#
#
# summaries = glob.glob("*.summary")
# n_seqs = list()
# for file in summaries:
#     # reading the file and discarding the "#" starting sequences returns this kind of line:
#     # 1     CL_1001                 39   338   275     0.80  0.589
#     lines = open(file).readlines()
#     for line in lines:
#         if not line.startswith("#") and not line.startswith("\n"):
#                 split = line.split(" ")
#                 while "" in split:
#                     split.remove("")
#
#                 # split is now like this:
#                 # ['1', 'CL_485', '30', '108', '90', '0.99', '0.635', '\n']
#                 n_seqs.append([f"{split[1]}.hmm", int(split[2])])
#
# order = sorted(n_seqs, key=lambda CL: int(CL[1]), reverse=True)
# print(order[:10])

##
hmms = sorted(glob.glob("*.hmm"))

# hmms_to_run = list()
#
# for CL in order:
#     # In case we want to run only .hmm with a certain number of sequences
#     #if CL[0] in hmms and CL[1] > 20:
#     if CL[0] in hmms:
#         hmms_to_run.append(CL[0])
#
# print(hmms_to_run[:10])

def run_hmm_hmmsearch(hmm_file):
    outdir = "../0_hmmsearch_raw"

    start = time.time()
    CL_id = os.path.basename(hmm_file).split(".hmm")[0]
    os.system(f"hmmsearch --cpu 5 --noali --notextw -o {outdir}/{CL_id}.hmmsearch --domtblout {outdir}/{CL_id}.domtblout {hmm_file} ../nr")
    print(CL_id)
    timeit(start)
    print()

print(f"{len(hmms)} profiles will be analyzed")


pool = multiprocessing.Pool(processes=14)
shared_content = pool.map(run_hmm_hmmsearch,hmms ) #hmms_to_run
pool.close()
pool.join()
