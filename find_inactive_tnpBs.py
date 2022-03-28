#!/usr/bin/env python
# coding: utf-8

# In[23]:


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


# In[3]:


#outdir = "../tnpBs/representatives"
outdir = "../tnpBs/all_tnpBs"


infile_tnpb_df = "tnpB_targetgenes_pfam.filtered.csv"
infile_tnpBs = "tnpBs_in_testdb.p100.1e4.txt"
tnpB_p60s_path = "../tnpBs/cluster/OUTPUT/tmp/clu_cluster.tsv"
infile_mafft = outdir + "/_all_tnpBs_mafft_input.faa"
outfile_mafft = outdir + "/_all_tnpBs_mafft.out.faa"

#outfile_tnpBs_inactive = "../tnpBs/_inactive_tnpBs.faa"
outfile_tnpBs_inactive_3 = "../tnpBs/_all_inactive_tnpBs.3.faa"
outfile_tnpBs_inactive_2 = "../tnpBs/_all_inactive_tnpBs.2.faa"


# In[25]:


tnpB_10_target_annot_df = pd.read_csv(infile_tnpb_df).iloc[:,1:]
#tnpB_10_target_annot_df = pd.read_csv("tnpB_targetgenes_pfam.csv").iloc[:,1:]


# In[26]:


tnpB_10_target_annot_df.shape


# In[11]:


tnpBs_list = []
with open(infile_tnpBs, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        p100 = line.split('\n')[0]
        tnpBs_list.append(p100)
tnpBs_set = set(tnpBs_list)
len(tnpBs_list), len(tnpBs_set)


# ## get high icity tnpBs (p60s) to align

# In[36]:


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
# # run only once to cluster all tnpBs in tnpBs_list
#get_faas_protidlist(tnpBs_list, "../tnpBs/cluster/INPUT/tnpB_mmseqs_input.faa")


# In[18]:


# run only after p60 clusters formed
tnpB_p60_df = pd.read_csv(tnpB_p60s_path, sep = '\t', header = None).rename(columns = {0:"p60",1:'p100'}).drop_duplicates()


# In[27]:


high_icity_tnpBs_ls = []
for baitp100s_str in tnpB_10_target_annot_df["baitp100s"]:
    baitp100s_ls = ast.literal_eval(baitp100s_str)
    high_icity_tnpBs_ls += baitp100s_ls
high_icity_tnpBs = [tnpB for tnpB in high_icity_tnpBs_ls if tnpB in tnpBs_set]
high_icity_tnpBs_set = set(high_icity_tnpBs)


# In[28]:


len(high_icity_tnpBs), len([tnpB for tnpB in high_icity_tnpBs_ls])


# In[29]:


tnpB_p60_list = list(set(tnpB_p60_df['p60']))
tnpB_p60_highicity_list = [p60 for p60 in tnpB_p60_list if p60 in high_icity_tnpBs_set]


# In[30]:


len(high_icity_tnpBs_set), len(tnpBs_list)


# In[32]:


len(tnpB_p60_highicity_list), len(tnpB_p60_list)


# ### run alignments and get catalytically inactive tnpBs

# In[37]:


def get_multifaa_protidlist(protidlist, outpath):
    with open(outpath, "w") as outfile:
        for protid in protidlist:
            protseq = get_prot_sequence(protid)
            if len(protseq) > 300 and len(protseq) < 500:
                print(">" + protid, file=outfile)
                print(protseq, file=outfile)
#get_multifaa_protidlist(tnpB_p60_highicity_list, infile_mafft)


# In[7]:


tnpBs_list = []
with open(infile_tnpBs, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        line=line.strip('\n')
        tnpBs_list.append(line)


# In[ ]:


# to do - delete later
#get_multifaa_protidlist(tnpBs_list, infile_mafft)


# In[40]:


cmd_mafft = "ginsi --thread 32 " + infile_mafft + " > " + outfile_mafft
print("running mafft ginsi: " + cmd_mafft)
os.system(cmd_mafft)


# In[42]:


ref_tnpB_pid = '0fea0aedf485f57c86'
seq1 = get_prot_sequence(ref_tnpB_pid)
seq1[187], seq1[271], seq1[354]


# In[47]:


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


# In[48]:





# In[99]:


count_3, count_2, count_1, count_0 = 0, 0, 0, 0
with open(outfile_mafft) as handle, open (outfile_tnpBs_inactive_3, "w") as outfile_3, open (outfile_tnpBs_inactive_2, "w") as outfile_2:
    for rec in SeqIO.parse(handle, 'fasta'):
        sequence = rec.seq
        pid = rec.id
        first_residue_missing = sequence[align_indices_dict[187]] != 'D'
        second_residue_missing = sequence[align_indices_dict[271]] != 'E'
        third_residue_missing = sequence[align_indices_dict[354]] != 'D'
        num_res_missing = first_residue_missing + second_residue_missing + third_residue_missing 
        if num_res_missing == 3:
            print(">" + pid, file = outfile)
            sequence_gapfree = get_prot_sequence(pid)
            print(sequence_gapfree, file = outfile)
            count_3 += 1
        elif num_res_missing == 2:
            print(">" + pid, file = outfile)
            sequence_gapfree = get_prot_sequence(pid)
            print(sequence_gapfree, file = outfile)
            count_2 += 1
        elif num_res_missing == 1:
            count_1 += 1
        elif num_res_missing == 0:
            count_0 += 1


# In[100]:


count_3, count_2, count_1, count_0


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
with open (outfile_tnpBs_inactive) as handle:
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




