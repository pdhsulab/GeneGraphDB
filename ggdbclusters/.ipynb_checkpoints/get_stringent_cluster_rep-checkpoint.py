#!/usr/bin/env python
# coding: utf-8

# In[42]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time


# In[117]:


# helper functions

#get protein viability from sql database. Input pid. Assume database conn is active
#con=sqlite3.connect("genegraph.db")
con=sqlite3.connect("80kprotein_stats.db")
cur = con.cursor()
def get_protein_is_goodrep(pid):
    #cmd = "SELECT * FROM proteins WHERE hashid='%s'" % pid
    cmd = "SELECT * FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return cur.fetchone()


# In[141]:


tic = time.time()
print(get_protein_is_goodrep("1854c240bfbd0005d449"))
toc = time.time()
#print(toc-tic)


# In[127]:




# In[163]:


# create function to get all the csv lines for a stringent cluster
# return value is list of lines to write to new file
def write_stringentcluster_outfile(pid_set, cur_rep, outfile):
    if len(pid_set) == 0:
        pass
    else:
        out_csv_row_list = [cur_rep + "," + pid for pid in pid_set]
        with open(outfile, "a") as f:
            for out_csv_row in out_csv_row_list:
                print(out_csv_row, file=f)
#get_stringent_rep(test_set, "3")
def return_better_rep(candidate_stats, cur_stats, cand_id, cur_id):
    cur_isgapfree, cur_length, cur_iscomplete = cur_stats[1], cur_stats[2], cur_stats[3]
    cand_isgapfree, cand_length, cand_iscomplete = candidate_stats[1], candidate_stats[2], candidate_stats[3]
    if not (cur_isgapfree == "True") & (cand_isgapfree == "True"):
        return cand_id
    elif (cur_isgapfree == "True") & (cand_isgapfree == "True"):
        if cand_length > cur_length:
            return cand_id
        else: return cur_id
    else:
        return cur_id
    
def log_badrepresentatives(rep_id, outfile):
    rep_stats = get_protein_is_goodrep(rep_id)
    rep_isgapfree, rep_iscomplete = (rep_stats[1] == "True"), (rep_stats[3] == "True")
    print(rep_isgapfree, rep_iscomplete)
    if (not rep_isgapfree) | (not rep_iscomplete):
        with open(outfile, "a") as f:
            print(rep_stats[1] + "," + rep_stats[3], file=f)


# In[130]:


# read in test clu_cluster.tsv file
#path_tsv = "../clusters/OUTPUT/stringent/tmp/clu_cluster_10k.tsv"
path_tsv = "../clusters/OUTPUT/stringent/tmp/clu_cluster.tsv"

#variables used to track shifting between clusters
prev_clusterid = ""
is_newgroup = True
tracker_ = 0
with open(path_tsv, "r") as in_file:
    for line in csv.reader(in_file, delimiter='\t'):
        if tracker_ == 1000:
            continue
        # track which group I am in
        cur_clusterid = line[0]
        cand_clusterrep = line[1][:20]
        is_newgroup = cur_clusterid != prev_clusterid
        if is_newgroup:
            if prev_clusterid != "": #if anything but first iteration
                pid_set.remove(cur_clusterrep)
                log_badrepresentatives(cur_clusterrep, "../clusters/clu_badrep_stringent.csv")
                write_stringentcluster_outfile(pid_set, cur_clusterrep, "../clusters/clu_rep_stringent_final.csv")
            else: 
                # create headers for log clu_badrep_stringent.csv
                with open("../clusters/clu_badrep_stringent.csv", "w") as f:
                    print("is_gap_free,iscomplete", file = f)
            cur_clusterrep = line[0][:20]
            cur_stats = get_protein_is_goodrep(cur_clusterrep)
            cand_stats = get_protein_is_goodrep(cand_clusterrep)
            pid_set = set()
            pid_set.add(cur_clusterrep)
            pid_set.add(cand_clusterrep)
            cur_clusterrep = return_better_rep(cand_stats, cur_stats, cand_clusterrep, cur_clusterrep)
            
        else:
            cur_stats = get_protein_is_goodrep(cur_clusterrep)
            cand_stats = get_protein_is_goodrep(cand_clusterrep)
            pid_set.add(cand_clusterrep)
            cur_clusterrep = return_better_rep(cand_stats, cur_stats, cand_clusterrep, cur_clusterrep) 
            ___
        
        # what if there is no representative protein that is complete? TO DO - record this to a csv. NOTE - fragments are already recorded in sql database. 
        
        # # write to new file
        # unique_pid = line[1]
        prev_clusterid = cur_clusterid
        tracker_ += 1
    


# In[134]:


# #data exploration of og clu_clusters file
# pd.read_csv(path_tsv, sep = "\t", header = None)
# unique_pid = pd.read_csv(path_tsv, sep = "\t", header = None)[1].unique()


# In[ ]:

con.close()


