#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv


# In[2]:


path_coords = "../clusters/INPUT/stringent/mmseqs2_testdb_input.faa"

outfile_stats = open("../clusters/proteins_full_stats.csv", "w")
print("pid,isgapfree,length,iscomplete,sequence", file=outfile_stats)
               
fasta_sequences = SeqIO.parse(open(path_coords),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    p_gapfree = sequence[:-1].count("*") == 0
    p_len = len(sequence)
    p_complete = (sequence[0] == "M")&(sequence[-1] == "*")
    print(name + "," + str(p_gapfree) + "," + str(p_len) + "," + str(p_complete) + "," + sequence, file=outfile_stats)


# In[ ]:




