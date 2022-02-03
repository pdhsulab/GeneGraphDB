#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv


# In[13]:


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


# In[14]:


def get_mmseqs_input(file_paths_list, outfile_path):
    num_samples = len(file_paths_list)
    final_cycle_ind = num_samples //1000
    cat_cmd_1 = "cat "
    for i in range(1000):
        cat_cmd_1 += file_paths_list[i] + " "
    cat_cmd_1 += " > " + outfile_path
    os.system(cat_cmd_1)

    for i in range(1, final_cycle_ind):
        cat_cmd = "cat "
        for j in range(i * 1000, (i + 1) * 1000):

            cat_cmd += file_paths_list[j] + " "
        cat_cmd += " >> " + outfile_path
        os.system(cat_cmd)

    cat_cmd = "cat "
    for i in range(final_cycle_ind * 1000, num_samples):
        cat_cmd += file_paths_list[i] + " "
    cat_cmd += " >> " + outfile_path
    os.system(cat_cmd)


# ### make multi prodigal .faa inputs

# In[15]:


prodigal_paths = []
for sample in drep_samples.keys():
    sample_path = drep_samples[sample] + sample
    prodigal_path = sample_path + "/" + sample + ".prodigal.faa.gz"
    prodigal_paths.append(prodigal_path)
    
outfile_path = "../clusters/INPUT/mmseqs2_testdb_input.faa.gz"

#get_mmseqs_input(prodigal_paths, outfile_path)


# In[16]:


backup_cmd = "gsutil cp " + outfile_path + " gs://jluo_bucket/ggdb/clusters"
#os.system(backup_cmd)


# ### make multi .gff inputs

# In[ ]:


print("making multi .gff inputs")
gff_paths = []
for sample in drep_samples.keys():
    sample_path = drep_samples[sample] + sample
    prodigal_path = sample_path + "/" + sample + ".prodigal.gff.gz"
    prodigal_paths.append(prodigal_path)
    
outfile_path_gff = "../clusters/mmseqs2_testdb_input.gff.gz"

get_mmseqs_input(prodigal_paths, outfile_path_gff)


# In[ ]:


backup_cmd_2 = "gsutil cp " + outfile_path_gff + " gs://jluo_bucket/ggdb/clusters"
#os.system(backup_cmd_2)


# ### data exploration of failed samples

# In[35]:


path_temp = "ggdb_multisql_errorlog.csv"


# In[ ]:





# In[47]:


in_df = pd.read_csv(path_temp, header=None)
in_df = in_df.iloc[3928:].set_index(0)
#in_df.to_csv(path_temp, header = False)


# ## choose representatives for stringent clusters

# In[ ]:





# In[7]:


# iterate through unique pid (cluster representatives) in tsv
# explore data (distribution of cluster sizes, do any representatives have >1 *)
# search for protein in fasta file (load in memory for now)
# length + # of *'s
# priority is (1) no * (2) length 

path_tsv = "../clusters/OUTPUT/stringent/tmp/clu_cluster_10k.tsv"
with open(path_tsv, "r") as in_file:
    for line in csv.reader(in_file, delimiter='\t'):
        # search for protein
        # get stats
        # write to new file
        unique_pid = line[1]

pid_to_sequence_dict = {}
path_coords = "../clusters/proteines_testdb_10k.faa"

# to do - this needs to be more scalable; can't load everything into memory
fasta_sequences = SeqIO.parse(open(path_coords),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    pid_to_sequence_dict[name] = sequence


def get_protein_info(pid):
    p_len = len(pid_to_sequence_dict[pid])


# ## Cluster_df data ex to get cluster summary stats

# In[22]:


cluster_df = pd.read_csv(path_tsv, sep = '\t', header = None)
rep_counts = cluster_df.groupby(0).count()[1]


# In[31]:


all_pids_s = cluster_df[1]


# In[ ]:




