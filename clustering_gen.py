#!/usr/bin/env python
# coding: utf-8

# In[10]:


import sys
import sqlite3
import subprocess
import os


# make mmseqs input multifaa, run mmseqs, get multifasta output, extract ids

# In[41]:


# pid_inpath, mmseqs_dir, identity = "cas1.txt", "../tnpBs/testclust", "0.6"
pid_inpath = sys.argv[1]
mmseqs_dir = sys.argv[2]
identity = sys.argv[3]
os.system("mkdir " + mmseqs_dir + " " + mmseqs_dir + "/INPUT " + mmseqs_dir + "/OUTPUT " + mmseqs_dir + "/OUTPUT/tmp")
mmseqs_faa_input = mmseqs_dir + "/INPUT/cluster_input.faa"
mmseqs_faa_output = mmseqs_dir + "/OUTPUT/"
mmseqs_faa_tmpoutput = mmseqs_dir + "/OUTPUT/tmp"


# ## get mmseqs input

# In[33]:


in_pid_list = []
with open(pid_inpath, "r") as infile:
    lines=infile.readlines()
    for line in lines:
        line = line.strip('\n')[:18]
        in_pid_list.append(line)


# In[35]:


con=sqlite3.connect("80kprotein_stats.db")
cur = con.cursor()
def get_prot_sequence(pid):
    cmd = "SELECT sequence FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    try:
        return str(cur.fetchone()[0])
    except:
        print(pid)
        asdf
def get_faas_protidlist(protidlist, outfile_path):
    with open(outfile_path, "w") as outfile:
        for protid in protidlist:
            protseq = get_prot_sequence(protid)
            print(">" + protid, file=outfile)
            print(protseq, file=outfile)
get_faas_protidlist(in_pid_list, mmseqs_faa_input)
con.close()


# In[42]:


cmd = "mmseqs easy-linclust " + mmseqs_faa_input + " " + mmseqs_faa_output + " " + mmseqs_faa_tmpoutput + " " + "--threads 32 -e 0.001 --min-seq-id " + identity + " -c " + identity + " --cov-mode 0 --spaced-kmer-mode 0 --remove-tmp-files 1"
process = subprocess.Popen(cmd.split(' '))
stdout, stderr = process.communicate()


# In[ ]:




