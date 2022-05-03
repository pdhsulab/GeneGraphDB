#!/usr/bin/env python
# coding: utf-8

# In[156]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SearchIO
from BCBio import GFF
import csv
import sqlite3
import time
from multiprocessing import Pool, cpu_count
import sys
from calc_icity_ import *
import ast
from collections import defaultdict
import subprocess
import numpy as np
import sys


# ### helper

# In[140]:


# make multi.faa file with all dTnpB p90 clusters
def get_prot_sequence(pid):
    con = sqlite3.connect("80kprotein_stats.db")
    cur = con.cursor()
    cmd = "SELECT sequence FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return str(cur.fetchone()[0]).replace("*", "")
    con.close()
def get_faas_protidlist(protidlist, outfile_path, minlen = 0):
    with open(outfile_path, "w") as outfile:
        for protid in protidlist:
            protseq = get_prot_sequence(protid)
            if len(protseq) > minlen:
                print(">" + protid, file=outfile)
                print(protseq, file=outfile)
def get_p90(bait_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p90 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)
def get_p30(bait_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p30 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)


# In[137]:


# import multi fasta file
inpath = sys.argv[1]
outpath = sys.argv[2]

input_name = inpath.split('/')[-1].replace(".faa", "")


# In[152]:


print(inpath)
print(outpath)


# ### kalign

# In[138]:


infile_k = inpath
outfile_k = "{}/{}_kalign.out.faa".format(outpath, input_name)

def kalign(infile, outfile):
    cmd = "kalign -i {} -o {}".format(infile, outfile)
    print(cmd)
    process = subprocess.Popen(cmd.split(' '))
    process.wait()
    


# In[ ]:


kalign(infile_k, outfile_k)


# ## Build tree

# In[ ]:


infile_tree = outfile_k
outfile_tree = outfile_k.replace("faa", "tree")
def fasttree(infile, outfile):
    cmd = "FastTree {} > {}".format(infile, outfile)
    print(cmd)
    time.sleep(2)
    os.system(cmd)
    #process = subprocess.Popen(cmd.split(' '))
fasttree(infile_tree, outfile_tree)


# ### scratch work - comment out later
# 

# #### make multi.faa file with all dTnpB p90 clusters

# In[109]:


with open("../ggdb_dfs/dTnpB_anymut_target_analysis.filtered.tsv", "r") as infile:
    df_ids = pd.read_csv(infile, sep = '\t').iloc[[0,1,2,14],:2]
df_ids["baitp90s"] = df_ids["baitp100s"].apply(lambda x: str([get_p90(pid) for pid in ast.literal_eval(x)]))
df_ids["baitp30s"] = df_ids["baitp100s"].apply(lambda x: str([get_p30(pid) for pid in ast.literal_eval(x)]))


# In[126]:


d_90, d_30 = defaultdict(set), defaultdict(set)


# In[141]:


for i in range(df_ids.shape[0]):
    target_id = df_ids.iloc[i,0]
    bait_p90ids = ast.literal_eval(df_ids.iloc[i,2])
    for pid in bait_p90ids:
        d_90[pid].update([target_id])
    bait_p30ids = ast.literal_eval(df_ids.iloc[i,3])
    for pid in bait_p30ids:
        d_30[pid].update([target_id])


# In[142]:


dTnpB_p90s = list(d_90.keys())
get_faas_protidlist(dTnpB_p90s, "../tree/sigma70/dTnpB_p90s.faa")


# In[143]:


dTnpB_p30s = list(d_30.keys())
get_faas_protidlist(dTnpB_p30s, "../tree/sigma70/dTnpB_p30s.faa")


# #### make multi.fna file with dTnpB locus
# 

# In[165]:


target_id = "23422406293d40c201"


# In[201]:


def get_dtnpb_faa(target_ids, outpath):
    for target_id in target_ids:
        gff_path = "../annotate/annotategff/OUTPUT/{}.sorted.gff".format(target_id)
        fna_path = "../annotate/annotategff/INPUT/fna/{}.examples.fna".format(target_id)
        with open(fna_path) as seq_handle:
            seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
        with open(gff_path) as in_handle, open(outpath, "w") as outfile:
            for rec in GFF.parse(in_handle):
                for i in range(len(rec.features)):
                    if rec.features[i].type == "CDS_bait":
                        bait_loc = rec.features[i].location
                        bait_id = rec.features[i].id
                        rec_id = rec.id
                        dTnpB_loci_seq = seq_dict[rec_id].seq[bait_loc.start + 0:bait_loc.end + 250]
                        rec_bait_id = "{}|{}".format(rec_id, bait_id)
                        print(">" + rec_bait_id, file=outfile)
                        print(dTnpB_loci_seq, file=outfile)


# In[200]:


target_ids = ["aa6af7e9289c3558d3", "23422406293d40c201", "c1050b21cc75640d51", "6bf6c4c7da68779d7a"]
outfile_dTnpB_fna = "{}/dTnpB_all_loci.fna".format(outpath)
get_dtnpb_faa(target_ids, outfile_dTnpB_fna)


# In[ ]:




