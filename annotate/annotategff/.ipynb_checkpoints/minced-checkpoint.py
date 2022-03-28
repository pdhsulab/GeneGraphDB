#!/usr/bin/env python
# coding: utf-8

# In[3]:


import subprocess
from multiprocessing import Pool, cpu_count


# In[2]:


def annot_crispr(tnpB_multifna_path, out_raw_path):
    cmd = "minced {} {}".format(tnpB_multifna_path, out_raw_path)
    process = subprocess.Popen(cmd.split(' '))


# In[ ]:


def annot_crispr_pool(ids_list):
    paths_list = []
    for targetid in ids_list:
        print("a")
        tnpB_multifna_path = 'INPUT/fna/' + targetid + '.examples.fna'
        print(targetid)
        out_raw_path = 'minced/' + targetid + '.minced.gff'
        paths_list.append((tnpB_multifna_path, out_raw_path))
    pool = Pool(cpu_count())
    results = pool.starmap(annot_crispr, iterable = paths_list)
    pool.close()
    pool.join()
    return results

