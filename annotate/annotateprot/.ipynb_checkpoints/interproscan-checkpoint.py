#!/usr/bin/env python
# coding: utf-8
# %%

# %%


from multiprocessing import Pool, cpu_count
from Bio import SearchIO
import os

def interproscan_pool(protein_ids):
    pool = Pool(cpu_count())
    results = pool.map(interproscan, iterable = protein_ids)
    pool.close()
    pool.join()
    return results

def interproscan(protein_id):
    infile_path = "INPUT/faa/" + protein_id + ".faa"
    outfile_path = "interproscan/" + protein_id + ".interproscan"
    cmd = "bash interproscan.sh -i " + infile_path + " -b " + outfile_path + " -f tsv"
    try:
        os.system(cmd)
    except:
        with open("interproscan/error_log.txt", "a") as outfile:
            print(protein_id, file=outfile)