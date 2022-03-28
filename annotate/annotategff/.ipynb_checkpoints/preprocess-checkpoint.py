#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import SeqIO
import os


# In[26]:


infile_paths = ['INPUT/gff/' + file for file in os.listdir('INPUT/gff') if file != '.ipynb_checkpoints']
infile_paths


# In[7]:


# integrate later
def single_gff_to_faa(infile_path):
    multifasta_path = infile_path.replace("gff", "faa")
    start_reading = False
    outfile_gff_path = infile_path.replace("INPUT/gff", "OUTPUT")
    with open(infile_path, 'r') as infile, open(outfile_gff_path, 'w') as outfile_gff, open(multifasta_path, "w") as outfile:
        for line in infile:
            line = line.strip('\n')
            if not start_reading:
                print(line, file=outfile_gff)
            if line.startswith(">"):
                start_reading = True
            if start_reading:
                print(line, file=outfile)


# In[ ]:




