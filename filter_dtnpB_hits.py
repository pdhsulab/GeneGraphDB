#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# ## helper functions
# def get_permissive_rep_(bait_pid):
# def get_p90s(p30_id):
# def get_p100s(p90_id):

# In[79]:


#### helper functions
def get_permissive_rep_(bait_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p30 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)

def get_p90s(p30_id):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd = "SELECT p90 FROM clusters WHERE p30 = '%s'" % (p30_id)
    cursor.execute(cmd)
    p90s = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(list(set(p90s)))

def get_p100s(p90_id):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd = "SELECT p100 FROM clusters WHERE p90 = '%s'" % (p90_id)
    cursor.execute(cmd)
    p100s = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(list(set(p100s)))


# In[40]:


tnpbdir = "../tnpBs/dtnpB_highicity"
inpath_dtnpB_targets = "dtnpB_target_fetchinput.csv"


# In[91]:


dtnpB_hiicity_ids = []
dtnpB_targets_df = pd.read_csv(inpath_dtnpB_targets)
dtnpB_series_p100 = dtnpB_targets_df["baitp100s"].apply(lambda x: ast.literal_eval(x))
for dtnpB_id_list in dtnpB_series_p100:
    dtnpB_hiicity_ids += dtnpB_id_list
dtnpB_hiicity_ids = list(set(dtnpB_hiicity_ids))
len(dtnpB_hiicity_ids)


# ## get high icity tnpBs (p60s) to align

# In[41]:


mmseqs_dir = tnpbdir
mmseqs_in = mmseqs_dir + "/in_ids"


# In[33]:


def cluster_(pid_inpath, mmseqs_dir, identity, pid_list = None):
    if pid_list != None:
        with open(pid_inpath, "w") as outfile:
            for pid in pid_list:
                print(pid, file=outfile)
    cmd = "python3 clustering_gen.py {} {} {}".format(pid_inpath, mmseqs_dir, identity)
    process = subprocess.Popen(cmd.split(' '))
    output_tsv = mmseqs_dir + "/OUTPUT/_cluster.tsv"
    cluster_rep_df = pd.read_csv(output_tsv, sep = '\t').drop_duplicates()
    return list(set(cluster_rep_df.iloc[:,0]))


# In[39]:


dtnpB_hiicity_p60_ids = cluster_(mmseqs_in, mmseqs_dir, "0.6", dtnpB_hiicity_ids)
len(dtnpB_hiicity_p60_ids)


# In[100]:


def p100_p60_cluster_dict(mmseqs_dir):
    cluster_dict = {}
    output_tsv = mmseqs_dir + "/OUTPUT/_cluster.tsv"
    cluster_rep_df = pd.read_csv(output_tsv, sep = '\t', header = None).drop_duplicates()
    for i in range(len(cluster_rep_df)):
        p60id = cluster_rep_df.iloc[i,0]
        p100id = cluster_rep_df.iloc[i,1]
        cluster_dict[p100id] = p60id
    return cluster_dict
p100_p60_dict = p100_p60_cluster_dict(mmseqs_dir)

def get_p60_rep(p100_list_str):
    p60_list = []
    p100_list = ast.literal_eval(p100_list_str)
    for p100 in p100_list:
        p60 = p100_p60_dict[p100]
        p60_list.append(p60)
    return p60_list


# In[105]:


dtnpB_targets_df["baitp60s"] = dtnpB_targets_df["baitp100s"].apply(lambda x: get_p60_rep(x))


# ### run alignments and get catalytically inactive tnpBs

# In[107]:


infile_mafft = tnpbdir + "/mafft_input.p60.faa"
outfile_mafft = tnpbdir + "/mafft_ouput.p60.faa"
infile_mafft_all = tnpbdir + "/mafft_input.p100.faa"
outfile_mafft_all = tnpbdir + "/mafft_ouput.p100.faa"


# In[37]:


def get_prot_sequence(pid):
    con=sqlite3.connect("80kprotein_stats.db")
    cur = con.cursor()
    cmd = "SELECT sequence FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return str(cur.fetchone()[0])
    con.close()


# In[108]:


def get_multifaa_protidlist(protidlist, outpath):
    with open(outpath, "w") as outfile:
        for protid in protidlist:
            protseq = get_prot_sequence(protid)
            if len(protseq) > 300 and len(protseq) < 500:
                print(">" + protid, file=outfile)
                print(protseq, file=outfile)
#get_multifaa_protidlist(dtnpB_hiicity_p60_ids + ["0fea0aedf485f57c86"], infile_mafft)
get_multifaa_protidlist(dtnpB_hiicity_ids + ["0fea0aedf485f57c86"], infile_mafft_all)


# In[48]:


# cmd_mafft = "ginsi --thread 32 " + infile_mafft + " > " + outfile_mafft
# print("running mafft ginsi: " + cmd_mafft)
# os.system(cmd_mafft)


# In[ ]:


cmd_mafft = "ginsi --thread 32 " + infile_mafft_all + " > " + outfile_mafft_all
print("running mafft ginsi: " + cmd_mafft)
os.system(cmd_mafft)
asdf


# In[49]:


ref_tnpB_pid = '0fea0aedf485f57c86'
ref_seq = get_prot_sequence(ref_tnpB_pid)
ref_seq[187], ref_seq[271], ref_seq[354]


# In[75]:


def get_alignindices_dict(ref_tnpB_pid, msa_file):
    # get seq1, which is the amino acid fasta format for ref_tnpB_pid
    seq1 = get_prot_sequence(ref_tnpB_pid)
    # get seq2, which is the amino acid fasta format with gaps, post-alignment
    with open(msa_file) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            pid = rec.id
            if pid == ref_tnpB_pid:
                seq2 = str(rec.seq)
                break
    # map the indices of amino acids between the two sequence formats
    index_to_alignindex = {}
    align_indexes = []
    for i in range(len(seq2)):
        char = seq2[i]
        if char != '-':
            align_indexes.append(i)
    for j in range(len(align_indexes)):
        index_to_alignindex[j] = align_indexes[j]
    return index_to_alignindex

align_indices_dict = get_alignindices_dict(ref_tnpB_pid, outfile_mafft)
align_indices_big_dict = get_alignindices_dict(ref_tnpB_pid, '../tmp.kalign.faa')


# In[ ]:


d1, e2, d3 = 0,0,0
with open(outfile_mafft) as handle:
    for rec in SeqIO.parse(handle, 'fasta'):
        p60_id = rec.id
        #if p60_id == ref_tnpB_pid:
        p60_seq = rec.seq
        residues = (p60_seq[align_indices_dict[187]], p60_seq[align_indices_dict[271]], p60_seq[align_indices_dict[354]])
        if residues[0] == "D":
            d1 += 1
        if residues[1] == "E":
            e2 += 1
        if residues[2] == "D":
            d3 += 1


# In[77]:


d1, e2, d3


# In[70]:


d1, e2, d3


# In[101]:


# p30_inactive_list = []
# with open(outfile_tnpBs_inactive) as handle:
#     for rec in SeqIO.parse(handle, 'fasta'):
#         p30_inactive = get_permissive_rep(rec.id)
#         p30_inactive_list.append(p30_inactive)
# p30_inactive_list.sort()
# set(p30_inactive_list)


# ### get hits near these inactive tnpBs

# In[103]:


inactive_tnpBs = []
with open(outfile_tnpBs_inactive) as handle:
    for rec in SeqIO.parse(handle, 'fasta'):
        sequence = rec.seq
        pid = rec.id
        inactive_tnpBs.append(pid)


# In[105]:


def extract_inactive_tnpBs(entry, inactive_tnpBs):
    baitp100s_set = set(ast.literal_eval(entry))
    intersection = baitp100s_set.intersection(inactive_tnpBs)
    return intersection


# In[106]:


def get_hits_with_inactive_tnpB_bait_df(inactive_tnpBs, df):
#     all_query_ids = [pfam for pfam in list(pfam_df['query_id'].sort_values().unique()) if type(pfam) == str]
#     query_ids = []
#     for query_ids_suffix in query_ids_suffixes:
#         query_ids += [pfam for pfam in all_query_ids if query_ids_suffix in pfam ]    
#     hits_s = df[df["query_id"].isin(query_ids)]['hit_id']
#     hits_set = set(hits_s)
#     hits_with_pfam_df = df[df["hit_id"].isin(hits_set)]
    
    
    baitlists_all = [ast.literal_eval(baitlist) for baitlist in df['baitp100s'].unique()]
    baitlists_expanded = []
    for baitlist in baitlists_all:
        for bait in baitlist:
            if bait in inactive_tnpBs:
                baitlists_expanded.append(str(baitlist))
                break
    hits_with_inactive_tnpB_df = df[df["baitp100s"].isin(baitlists_expanded)].dropna()
    # add column, "inactive tnpBs"
    hits_with_inactive_tnpB_df['inactive tnpBs'] = hits_with_inactive_tnpB_df['baitp100s'].apply(lambda x: extract_inactive_tnpBs(x, inactive_tnpBs))
    tnpB_hits_coloc_inactive_tnpB =  hits_with_inactive_tnpB_df[hits_with_inactive_tnpB_df['is_tnpB'] == 1]
                
    print("number hits with >1 inactive tnpB bait: " + str(len(hits_with_inactive_tnpB_df['hit_id'].unique())))
    print("number tnpB hits with >1 inactive tnpB bait: " + str(len(tnpB_hits_coloc_inactive_tnpB['hit_id'].unique())))
    return hits_with_inactive_tnpB_df


# In[107]:


pd.set_option('display.max_rows', None)
targets_inactivebaits_df = get_hits_with_inactive_tnpB_bait_df(inactive_tnpBs, tnpB_10_target_annot_df)
targets_inactivebaits_df["hit_id"].unique()


# In[108]:


target_to_bait_df = targets_inactivebaits_df[['hit_id', 'inactive tnpBs']]


# In[109]:


target_to_bait_df[target_to_bait_df["hit_id"] == '23422406293d40c201'].iloc[1, 1]


# In[115]:


'7f51222e4a7ac0fcf6' in high_icity_tnpBs_set


# In[ ]:




