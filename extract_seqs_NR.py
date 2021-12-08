'''
After parsing the hmmsearch results, extracts the sequences from the NR .faa file
'''

from Bio import SeqIO, SearchIO
import os, glob, argparse, multiprocessing, time
from functools import partial



def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-i', '--hmmsearch_dir',
                               dest='hmmsearch_dir',
                               required=True,
                               help='directory with the hmmsearch results'
                               )
    requiredArgs.add_argument('-f', '--nr_faa',
                               dest='nr_faa',
                               required=True,
                               help='NR.faa file'
                               )
    requiredArgs.add_argument('-o', '--outdir',
                               dest='outdir',
                               required=True,
                               help='directory to save the .faa files'
                               )

    return parser.parse_args()

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

def parse_hmmsearch_results(hmmsearch_file):
    cl_id = os.path.basename(hmmsearch_file).split(".")[0]
    #print(cl_id, os.path.getsize(hmmsearch_file))
    hits = [cl_id]
    # with open(hmmsearch_file, "r") as handle:
    #     for record in SearchIO.parse(handle, "hmmer3-text"):
    #         for hit in record.hits:
    #             for hsp in hit.hsps:
    #                 if hsp.evalue < 0.001:
    #                     hits.append(hit.id)
    #                     break

    # I write this to try to speed up the filtering process by removing all the hits with evalue > 0.001
    records = SearchIO.parse(hmmsearch_file, "hmmer3-text")
    all_hits = [hit for record in records for hit in record.hits if hit.evalue < 0.001]
    for hit in all_hits:
        #print(cl_id, hit.id)
        for hsp in hit.hsps:
            if hsp.evalue < 0.001:
                hits.append(hit.id)
                break

    return hits




def main():

    args = parse_args()

    # list the hmmsearch results files
    hmmsearch_files = glob.glob(f"{args.hmmsearch_dir}/*.hmmsearch")

    # remove already processed hmmsearch files
    done_cls = [os.path.basename(file).split(".")[0] for file in glob.glob(f"{args.outdir}/*.faa")]
    print(f"{len(done_cls)} already extracted sequences.")

    exclude = ["cl_s_321"]
    hmmsearch_files = [hmmsearch_file for hmmsearch_file in hmmsearch_files
                       if os.path.basename(hmmsearch_file).split(".")[0] not in done_cls + exclude]

    print(f"{len(hmmsearch_files)} families will be processed to extract hits sequences")


    print("parsing hmmsearch files...")
    start = time.time()
    # process hmmsearch files in parallel
    # cls_seqsids is a list of lists, in each list the first item is the cl_id, the rest are hits protein ids
    pool = multiprocessing.Pool(processes=70)
    cls_seqsids = pool.map(parse_hmmsearch_results, hmmsearch_files)
    pool.close()
    pool.join()
    print("Done")
    timeit(start)
    print()

    print("Indexing NR faa file...")
    start = time.time()
    nr_records = SeqIO.index(args.nr_faa, "fasta")
    print("Done")
    timeit(start)
    print()

    print(f"Extracting sequences for {len(cls_seqsids)} CLs...")
    start = time.time()
    for hit_ids in cls_seqsids:
        cl_id = hit_ids[0]
        #print(cl_id)
        outfile = f"{args.outdir}/{cl_id}.faa"
        hit_ids.pop(0)
        if hit_ids:
            #print(hit_ids)
            to_write = [nr_records[hit_id] for hit_id in hit_ids]
            with open(outfile, "w") as fout:
                SeqIO.write(to_write, fout, "fasta")
        else:
            os.system(f"touch {outfile}")

    print("Done")
    timeit(start)

if __name__ == "__main__":
    main()
