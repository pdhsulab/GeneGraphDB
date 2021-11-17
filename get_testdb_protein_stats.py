#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3


# In[10]:


# create proteins_full_stats.csv


# In[2]:


path_coords = "../clusters/INPUT/stringent/mmseqs2_testdb_input.faa"

outfile_stats = open("../clusters/proteins_full_stats.csv", "w")
print("pid,isgapfree,length,iscomplete,sequence", file=outfile_stats)
               
fasta_sequences = SeqIO.parse(open(path_coords),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id[:20], str(fasta.seq)
    p_gapfree = sequence[:-1].count("*") == 0
    p_len = len(sequence)
    p_complete = (sequence[0] == "M")&(sequence[-1] == "*")
    print(name + "," + str(p_gapfree) + "," + str(p_len) + "," + str(p_complete) + "," + sequence, file=outfile_stats)


# In[3]:


#load into sql database


# In[9]:


con = sqlite3.connect('80kprotein_stats.db')
cur = con.cursor()
cur.execute('''CREATE TABLE proteins (pid text, isgapfree integer, length real, iscomplete integer, sequence text, PRIMARY KEY(pid))''')
protein_csv = open(path_coords)
protein_csv = open("../clusters/proteins_full_stats.csv")
rows = csv.reader(protein_csv)
next(rows)
cmd = '''
INSERT OR IGNORE INTO proteins (pid, isgapfree, length, iscomplete, sequence) VALUES (?,?,?,?,?)
'''
cur.executemany(cmd, rows)
con.commit()
con.close()

