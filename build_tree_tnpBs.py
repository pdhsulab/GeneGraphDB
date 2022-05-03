#!/usr/bin/env python
# coding: utf-8

# In[45]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SearchIO
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


# In[83]:


def get_prot_sequence(pid):
    con=sqlite3.connect("80kprotein_stats.db")
    cur = con.cursor()
    cmd = "SELECT sequence FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return str(cur.fetchone()[0])
    con.close()
def get_faas_protidlist(protidlist, outfile_path):
    with open(outfile_path, "w") as outfile:
        for protid in protidlist:
            protseq = get_prot_sequence(protid)
            print(">" + protid, file=outfile)
            print(protseq, file=outfile)
    


# ### cluster

# In[ ]:


inactive_tnpBs_list = []
with open("../inactive_tnpbs_ref.3.faa", "r") as infile:
    lines=infile.readlines()
    for line in lines:
        if line[0] == '>':
            tnpBid = line.strip('>').strip('\n')
            inactive_tnpBs_list.append(tnpBid)


# In[ ]:


pid_inpath = "../tnpBs/_all_inactive_tnpBs.3.txt"
mmseqs_dir = "../tnpBs/clust_inactive"
identity = "0.6"
pid_list = inactive_tnpBs_list


# In[ ]:


def cluster_(pid_inpath, mmseqs_dir, identity, pid_list = None):
    if pid_list != None:
        with open(pid_inpath, "w") as outfile:
            for pid in pid_list:
                print(pid, file=outfile)
    cmd = "python3 clustering_gen.py {} {} {}".format(pid_inpath, mmseqs_dir, identity)
    process = subprocess.Popen(cmd.split(' '))
    output_tsv = mmseqs_dir + "/OUTPUT/_cluster.tsv"
    cluster_rep_df = pd.read_csv(output_tsv, sep = '\t')
    return list(set(cluster_rep_df.iloc[:,0]))


# In[ ]:


dTnpB_p60_list = cluster_(pid_inpath, mmseqs_dir, identity, pid_list = None)


# In[ ]:


len(pid_list), len(dTnpB_p60_list)


# ### make alignment multifaa infile

# - always start with all p60 tnpBs id.txt
# - add 100 high icity dtnpBs to "ignore" set
# - add 100 low icity dtnpBs to "ignore" set
# - select 300 random p60 tnpB ids

# In[88]:


# # concatenate catalytically active tnpBs to inactive list + condense for visibility
tnpB_p60s_path = "../tnpBs/cluster/OUTPUT/tmp/clu_cluster.tsv"
tnpB_p60_df = pd.read_csv(tnpB_p60s_path, sep = '\t', header = None).rename(columns = {0:"p60",1:'p100'}).drop_duplicates()
tnpB_p60_set = set(tnpB_p60_df["p60"])

tnpB_df = pd.read_csv("tnpB_icity_output.tsv", sep='\t')
tnpB_df_low_icity = tnpB_df[tnpB_df["denom"] > 10][tnpB_df["icity"] < .7].sort_values(["icity","numer"], ascending = False).drop_duplicates()
tnpB_df = tnpB_df[tnpB_df["icity"] > .7].sort_values(["icity","numer"], ascending = False).drop_duplicates()
tnpB_10_target_annot_df = tnpB_df[tnpB_df["denom"] > 10]

high_icity_dtnpBs_ls = []
for baitp100s_str in tnpB_10_target_annot_df["baitp100s"]:
    baitp100s_ls = ast.literal_eval(baitp100s_str)
    high_icity_dtnpBs_ls += baitp100s_ls
high_icity_dtnpBs_ls = list(set(high_icity_dtnpBs_ls))
high_icity_dtnpBs_samp = set(np.random.choice(high_icity_dtnpBs_ls, 100)).intersection(tnpB_p60_set)

low_icity_dtnpBs_ls = []
for baitp100s_str in tnpB_df_low_icity["baitp100s"]:
    baitp100s_ls = ast.literal_eval(baitp100s_str)
    low_icity_dtnpBs_ls += baitp100s_ls
low_icity_dtnpBs_ls = list(set(low_icity_dtnpBs_ls))
low_icity_dtnpBs_samp = set(np.random.choice(low_icity_dtnpBs_ls, 100))

tnpBs_p60samp = [tnpB for tnpB in np.random.choice(
    list(set(tnpB_p60_df["p60"])), 300, replace = False) 
                 if (tnpB not in low_icity_dtnpBs_samp) 
                 and (tnpB not in high_icity_dtnpBs_samp)]


# In[104]:


pid_inpath = "../tnpBs/_all_inactive_tnpBs.3.txt"
mmseqs_dir = "../tnpBs/clust_inactive"
identity = "0.6"
dtnpBs_samp_set = high_icity_dtnpBs_samp.update(low_icity_dtnpBs_samp)
#pid_list = list(dtnpBs_samp_set.update(set(tnpBs_p60samp)))
dtnpBs_samp_set


# In[116]:


high_icity_dtnpBs_samp.update({2})


# In[87]:


protidlist = 
kalign_infile = "../tnpBs/alignments/kalign_dTnpB_p60_ref.in.faa"
get_faas_protidlist(protidlist, kalign_infile)


# In[78]:


low_icity_dtnpBs_samp.intersection(high_icity_dtnpBs_samp)


# ### kalign

# In[25]:


infile_k = kalign_infile
outfile_k = "../tnpBs/alignments/kalign_dTnpB_p60_ref.out.faa"


# In[26]:


def kalign(infile, outfile, idlist = None):
    if idlist != None:
        with open(infile, "w") as kalign_infile:
            for protid in idlist:
                protseq = get_prot_sequence(protid)
                print(">" + protid, file=kalign_infile)
                print(protseq, file=kalign_infile)
    cmd = "kalign -i {} -o {}".format(infile, outfile)
    process = subprocess.Popen(cmd.split(' '))


# In[ ]:


kalign(infile_k, outfile_k, inactive_tnpBs_list)


# ### fasttree

# In[29]:


infile_tree = outfile_k
outfile_tree = outfile_k.replace("faa", "tree").replace("kalign_", "")


# In[30]:


def fasttree(infile, outfile):
    cmd = "FastTree {} > {}".format(infile, outfile)
    print(cmd)
    process = subprocess.Popen(cmd.split(' '))


# In[24]:


#fasttree(infile_tree, outfile_tree)


# In[ ]:




