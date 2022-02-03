#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time


# In[19]:


stats_outfile_name = "cluster_stats.txt"
stats_outfile = open(stats_outfile_name, "a")
# os.system("touch %s" % stats_outfile)


# In[ ]:


con=sqlite3.connect("80kprotein_stats.db")
cur = con.cursor()
def get_protein_is_goodrep(pid):
    #cmd = "SELECT * FROM proteins WHERE hashid='%s'" % pid
    cmd = "SELECT * FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return str(cur.fetchone())


# In[18]:


# with open(path_stringent_repseq) as f:
#     pass


# In[17]:


path_stringent_repseq = "../clusters/OUTPUT/stringent/tmp/clu_rep_seq.fasta"
path_permissive_repseq = "../clusters/OUTPUT/permissive/tmp/clu_rep_seq.fasta"


# In[ ]:


# get number clusters
def get_number_entries(fasta_path):
    cmd = "grep -c '>' " + fasta_path
    rv = os.popen(cmd).read()
    return rv

num_stringent = get_number_entries(path_stringent_repseq)
num_perm = get_number_entries(path_permissive_repseq)
print(num_stringent, num_perm)
print("number stringent clusters is " + num_stringent, file = stats_outfile)
print("number permissive clusters is " + num_perm, file = stats_outfile)


# In[20]:


# get distribution of genes per cluster
#def get_genespercluster_dist():
    


# In[ ]:


# skeleton code to save figures
# fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
# ax.plot([0,1,2], [10,20,3])
# fig.savefig('path/to/save/image/to.png')   # save the figure to file
# plt.close(fig)


# In[ ]:


con.commit()
con.close()

