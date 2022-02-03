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


# In[ ]:





# In[3]:


def get_permissive_rep(stringent_rep):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd_p = "SELECT reppid FROM permissive WHERE pid = '%s'" % (stringent_rep)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)


# In[4]:


def get_protein_seq(pid):
    con=sqlite3.connect("80kprotein_stats.db")
    cur = con.cursor()
    #cmd = "SELECT * FROM proteins WHERE hashid='%s'" % pid
    cmd = "SELECT * FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return cur.fetchone()[-1]
    con.close()


# In[ ]:


path_clu_rep = "../clusters/clu_rep_stringent_final.csv"


# In[25]:


outfile = open("../clusters/INPUT/permissive/clu_perm_mmseqs_input.faa", "w")
with open(path_clu_rep, 'r') as file:
    reader = csv.reader(file)
    next(reader)
    prev_rep = ""
    for row in reader:
        stringent_rep = row[0][:18]
        if prev_rep != stringent_rep:
            print(">" + stringent_rep, file = outfile)
            print(get_protein_seq(stringent_rep), file = outfile)
        prev_rep = stringent_rep


# In[ ]:




