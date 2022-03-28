#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
from multiprocessing import Pool


# In[22]:


infile_paths = ['INPUT/faa/' + file for file in os.listdir('INPUT/faa') if file != '.ipynb_checkpoints']
p30s = [infile_path.split('/')[-1].split('.')[0] for infile_path in infile_paths]


# In[24]:


p30s, infile_paths


# In[27]:


def cctyper(infile_paths):
    os.system("bash activate_cct.sh")
    for infile_path in infile_paths:
        p30 = infile_path.split('/')[-1].split('.')[0]
        outfile_path = 'cct/OUTPUT/' + p30
        cmd = "cctyper " + infile_path + " " + outfile_path + " --prodigal meta"
        os.system(cmd)
    os.system("bash activate_ggdb.sh")
cctyper(infile_paths)


# In[ ]:




