#!/usr/bin/env python
# coding: utf-8

# In[11]:


import os
import preprocess
import pfamannot
import parse


# ### for each .examples.gff file (target p30 neighborhood)
# - get multifasta file from gff: input for hmmsearch and cctyper
# - parse hmmsearch output
# - parse cctyper outputs
# 
# ### repeat for each .examples.gff file in input directory
# - 
# 

# In[6]:


# load input files into input directory
hits_near_inactive_tnpBs = ['7205dff302ff900300', '6c3b2ade31d3833745', '23422406293d40c201',
       'c1050b21cc75640d51', 'aa6af7e9289c3558d3', 'd2246f26c16fb9eec9',
       '6bf6c4c7da68779d7a']
input_files = [hitid + ".examples.gff" for hitid in hits_near_inactive_tnpBs]
input_filepaths = ["../ggdbfetch_output/" + hitpath for hitpath in input_files]
for i in range(len(input_filepaths)):
    file_path = input_filepaths[i]
    cmd = "cp " + file_path + " INPUT/gff/" + input_files[i]
    os.system(cmd)


# In[7]:


infile_paths = ['INPUT/gff/' + file for file in os.listdir('INPUT/gff') if file != '.ipynb_checkpoints']


# In[8]:


infile_paths


# In[ ]:





# In[14]:


def update_gffs(ids_list):
    create_faas(infile_paths)
    #annotate()
    #pfamannot.hmmsearch_pool(ids_list)
    final_gff_output(ids_list)
    
def create_faas(infile_paths):
    for path in infile_paths:
        preprocess.single_gff_to_faa(path)

def final_gff_output(ids_list):
    for hitid in ids_list:
        parse.put_gff_together(hitid)
        
def annotate():
    os.system('conda run -n cctyper python cctyper.py')
    


# In[ ]:


#annotate()


# In[13]:


update_gffs(hits_near_inactive_tnpBs)


# In[ ]:




