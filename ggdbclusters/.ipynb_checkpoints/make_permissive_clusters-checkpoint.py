#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time


# In[35]:


con=sqlite3.connect("80kprotein_stats.db")
cur = con.cursor()
def get_protein_seq(pid):
    #cmd = "SELECT * FROM proteins WHERE hashid='%s'" % pid
    cmd = "SELECT * FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return cur.fetchone()[-1]


# In[4]:


path_clu_rep = "../clusters/clu_rep_stringent_final.csv"


# In[56]:


outfile = open("../clusters/INPUT/permissive/clu_perm_mmseqs_input.faa", "w")
with open(path_clu_rep, 'r') as file:
    reader = csv.reader(file)
    next(reader)
    prev_rep = ""
    for row in reader:
        stringent_rep = row[0]
        if prev_rep != stringent_rep:
            print(">" + stringent_rep, file = outfile)
            print(get_protein_seq(stringent_rep), file = outfile)
        prev_rep = stringent_rep


# In[ ]:


con.commit()
con.close()

