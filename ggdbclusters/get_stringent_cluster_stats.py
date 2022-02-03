#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time


# In[6]:


stats_outfile_name = "cluster_stats.txt"
# os.system("touch %s" % stats_outfile_name)


# In[7]:


def get_protein_is_goodrep(pid):
    con=sqlite3.connect("80kprotein_stats.db")
    cur = con.cursor()
    cmd = "SELECT * FROM proteins WHERE pid = '%s'" % pid 
    cur.execute(cmd)
    return str(cur.fetchone())
    con.close()


# In[9]:


path_stringent_repseq = "../clusters/INPUT/permissive/clu_perm_mmseqs_input.faa"
path_permissive_repseq = "../clusters/OUTPUT/permissive/tmp/clu_rep_seq.fasta"


# In[7]:


# # get number clusters
# def get_number_entries(fasta_path):
#     cmd = "grep -c '>' " + fasta_path
#     rv = os.popen(cmd).read()
#     return rv
# print("number stringent clusters is " + get_number_entries(path_stringent_repseq), file = stats_outfile)
# print("number permissive clusters is " + get_number_entries(path_permissive_repseq), file = stats_outfile)


# In[10]:


# # get distribution of genes per cluster
# def get_genespercluster_dist(tsv_path):
#     df = pd.read_csv(tsv_path, sep="\t", header = None)
#     df = df.groupby(0).count()
#     fig, ax = plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
#     (n, bins, patches) = ax.hist(df[1], bins = [1,2,3,4,5,6,7,8,10,14,20,40,100,200,500,1000])
#     return str(n) + str(bins)
path_stringent_tsv = "../clusters/clu_rep_stringent_final.csv"
path_permissive_tsv = "../clusters/OUTPUT/permissive/tmp/clu_cluster.tsv"
# with open(stats_outfile_name, "a") as stats_outfile:
#     print("stats on genes per stringent cluster: " + get_genespercluster_dist(path_stringent_tsv), file = stats_outfile)
#     print("stats on genes per permissive cluster: " + get_genespercluster_dist(path_permissive_tsv), file = stats_outfile)


# In[ ]:


# get number fragments overall
# def get_num_fragments():
#     cmd = "SELECT count(*) FROM proteins WHERE iscomplete = 'False'"
#     cur.execute(cmd)
#     return str(cur.fetchone())
# with open(stats_outfile_name, "a") as stats_outfile:
#     num_frag = get_num_fragments()
#     print(num_frag)
#     print("number of fragments is: " + num_frag, file = stats_outfile)


# ### create 

# In[ ]:




# get number of singleton permissive/stringent clusters represented by fragments
def create_sqldb_fromclustertsv(tsv_path, is_perm):
    conn = sqlite3.connect('clusters.db')
    cursor = conn.cursor()
    if is_perm:
        cursor.execute('''CREATE TABLE permissive (reppid text, pid text, PRIMARY KEY(pid))''')
        protein_tsv = open(tsv_path, "r")
        rows = csv.reader(protein_tsv, delimiter="\t")
        cmd = '''
        INSERT OR IGNORE INTO permissive (reppid, pid) VALUES (?,?)
        '''
    else:
        cursor.execute('''CREATE TABLE stringent (reppid text, pid text, PRIMARY KEY(pid))''')
        protein_tsv = open(tsv_path, "r")
        rows = csv.reader(protein_tsv)
        cmd = '''
        INSERT OR IGNORE INTO stringent (reppid, pid) VALUES (?,?)
        '''
    cursor.executemany(cmd, rows)
    conn.commit()
    conn.close()

# create_sqldb_fromclustertsv(path_permissive_tsv, True)

#create_sqldb_fromclustertsv(path_stringent_tsv, False)

# def merge_80kproteinstats_permissive():
#     conn = sqlite3.connect('clusters.db')
#     cursor = conn.cursor()
#     #cmd1 = "attach 'clusters.db' as db1"
#     cmd2 = "ATTACH DATABASE '80kprotein_stats.db' AS db2"
#     cmd3 = "CREATE TABLE permissive_stats (reppid text, pid text, isgapfree integer, length real, iscomplete integer, sequence text, PRIMARY KEY(pid))"
#     cmd4 = "INSERT INTO permissive_stats SELECT reppid, a.pid, isgapfree, length, iscomplete, sequence FROM permissive a INNER JOIN db2.proteins b ON a.pid = b.pid"
#     #cursor.execute(cmd1)
#     conn.execute(cmd2)
#     #cursor.execute(cmd3)
#     #cursor.execute("SELECT * FROM db2.proteins LIMIT 10")
#     cursor.execute(cmd4)
#     conn.commit()
#     conn.close()
    
# merge_80kproteinstats_permissive()

def merge_80kproteinstats_stringent():
    conn = sqlite3.connect('clusters.db')
    cursor = conn.cursor()
    #cmd1 = "attach 'clusters.db' as db1"
    cmd2 = "ATTACH DATABASE '80kprotein_stats.db' AS db2"
    cmd3 = "CREATE TABLE stringent_stats (reppid text, pid text, isgapfree integer, length real, iscomplete integer, sequence text, PRIMARY KEY(pid))"
    cmd4 = "INSERT INTO stringent_stats SELECT reppid, a.pid, isgapfree, length, iscomplete, sequence FROM stringent_ a INNER JOIN db2.proteins b ON a.pid = b.pid"
    #cursor.execute(cmd1)
    conn.execute(cmd2)
    cursor.execute(cmd3)
    #cursor.execute("SELECT * FROM db2.proteins LIMIT 10")
    cursor.execute(cmd4)
    conn.commit()
    conn.close()

    
## Uncomment this
merge_80kproteinstats_stringent()


# In[22]:


# test_df = pd.read_csv("../clusters/stringent_stats.csv", nrows=200)
# test_df.to_csv("test.csv")
# test_df = test_df.set_index("reppid")
# test_df["temp"] = test_df.groupby("reppid")['length'].apply(lambda x: (max(x)-min(x))/max(x))
# test_df.reset_index()


# In[29]:


# # Final stats that require protein clusters merged with protein stats
# def get_length_distance_stats(path_csv):
#     df = pd.read_csv(path_csv)
#     df = df.set_index("reppid")
#     df["len_coverage"] = df.groupby("reppid")['length'].apply(lambda x: (max(x)-min(x))/max(x))
#     df = df.reset_index()
#     df[df["len_coverage"] > .2].to_csv("../clusters/clusters_poorcoverage.csv")
#     fig, ax = plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
#     (n, bins, patches) = ax.hist(df["len_coverage"], bins = [.02,.04,.06,.08,.1,.12,.14,.16,.18,.2,.3,1])

#     with open(stats_outfile_name, "a") as stats_outfile:
#         print("distribution of length coverage of fragments across clusters: " + str(n), file = stats_outfile)
#         print("distributed across bins: " + str(bins), file = stats_outfile)
        
# # get_length_distance_stats("../clusters/stringent_stats.csv")
# get_length_distance_stats("../clusters/permissive_stats.csv")

# def count_singleton_p_clusters_fragments():
#     conn = sqlite3.connect('clusters.db')
#     cursor = conn.cursor()
#     cmd_perm = "SELECT count(*) FROM (SELECT * FROM permissive_stats WHERE iscomplete = 'False' GROUP BY reppid HAVING count(reppid) = 1)"
#     cursor.execute(cmd_perm)
#     num_frag_perm = str(cursor.fetchone())
#     print(num_frag_perm)
#     with open(stats_outfile_name, "a") as stats_outfile:
#         print("number of perm fragments in singleton clusters is: " + num_frag_perm, file = stats_outfile)
# def count_singleton_s_clusters_fragments():
#     conn = sqlite3.connect('clusters.db')
#     cursor = conn.cursor()
#     cmd_stringent = "SELECT count(*) FROM (SELECT * FROM stringent_stats WHERE iscomplete = 'False' GROUP BY reppid HAVING count(reppid) = 1)"
#     cursor.execute(cmd_stringent)
#     num_frag_stringent = str(cursor.fetchone())
#     print(num_frag_stringent)
#     with open(stats_outfile_name, "a") as stats_outfile:
#         print("number of stringent fragments in singleton clusters is: " + num_frag_stringent, file = stats_outfile)
# count_singleton_p_clusters_fragments()
# count_singleton_s_clusters_fragments()


# In[ ]:


# skeleton code to save figures
# fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
# ax.plot([0,1,2], [10,20,3])
# fig.savefig('path/to/save/image/to.png')   # save the figure to file
# plt.close(fig)


# In[ ]:




