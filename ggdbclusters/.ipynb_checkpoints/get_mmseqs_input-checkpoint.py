#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os


# In[3]:


path_drep = "../drep_genomes/OUTPUT/rep_genomes/"
drep_samples = {}
for directory in os.listdir(path_drep):
    path_2 = path_drep + directory
    if os.path.isdir(path_2):
        for directory_2 in os.listdir(path_2):
            path_3 = path_2 + "/" + directory_2
            if os.path.isdir(path_3):
                for directory_3 in os.listdir(path_3):
                    path_4 = path_3 + "/" + directory_3 + "/"
                    if os.path.isdir(path_4):
                        for directory_samp in os.listdir(path_4):
                            samp_dir = path_4 + directory_samp
                            if os.path.isdir(samp_dir):
                                drep_samples[directory_samp] = path_4


# In[4]:


prodigal_paths = []
for sample in drep_samples.keys():
    sample_path = drep_samples[sample] + sample
    prodigal_path = sample_path + "/" + sample + ".prodigal.faa.gz"
    prodigal_paths.append(prodigal_path)


# In[ ]:


path_input = "../clusters/INPUT/mmseqs2_testdb_input.faa.gz"
os.system("cat " + prodigal_paths[0] + " > " + path_input)


# In[11]:


cat_cmd = "cat "
for i in range(1, len(prodigal_paths)):
    cat_cmd += prodigal_paths[i] + " "
cat_cmd += " >> " + path_input


# In[ ]:


os.system(cat_cmd)


# In[13]:


backup_cmd = "gsutil cp " + path_input + " gs://jluo_bucket/ggdb/clusters"
os.system(backup_cmd)


# In[ ]:




