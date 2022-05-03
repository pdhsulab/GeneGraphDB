#!/usr/bin/env python
# coding: utf-8

# In[26]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time
from multiprocessing import Pool, cpu_count
import sys
import ast
from collections import defaultdict


# In[ ]:





# ### simple -icity calculation: given single protein, find bait neighbourhoods of all related proteins within the 5kb window

# In[27]:


def get_list_ids_fromcursor(fetchall):
    return [fetchone[0] for fetchone in fetchall]


# In[51]:


def get_permissive_rep_(bait_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p30 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)


# In[28]:


def get_permissive_rep(bait_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p30 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    try:
        perm_rep = cursor.fetchone()[0]
    except:
        pass
    #perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)

def get_related_baits(perm_rep):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd_get_strin_rep = "SELECT p100 FROM clusters WHERE p30 = '%s'" % (perm_rep)
    cursor.execute(cmd_get_strin_rep)
    p100s_ls = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(list(set(p100s_ls)))


# In[29]:


def get_bait_neighbourhood(in_pid):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()

    neighbours_set = set()
    cmd_getbaitneighbs = "SELECT p2hash FROM prot2protwindow WHERE p1hash = '%s'" % (in_pid)
    #print(cmd_getbaitneighbs)
    cursor.execute(cmd_getbaitneighbs)
    for neighb_id in get_list_ids_fromcursor(cursor.fetchall()):
        neighbours_set.add(neighb_id) 
    cmd_getbaitneighbs = "SELECT p1hash FROM prot2protwindow WHERE p2hash = '%s'" % (in_pid)
    cursor.execute(cmd_getbaitneighbs)
    for neighb_id in get_list_ids_fromcursor(cursor.fetchall()):
        neighbours_set.add(neighb_id)
    conn.close()
    rv = defaultdict(list)
    for neighborid in neighbours_set:
        rv[neighborid].append(in_pid)
    return rv

# #this function was only used for calculating cas1/cas2 icity as positive controls
def get_full_bait_neighbourhood(baits_list):
    neighbours_set = set()
    for bait_id in baits_list:
        neighbours_set.update(get_bait_neighbourhood(bait_id))
    with open("bait_neighbourhood.txt", "w") as outfile:
        for neighbour_id in neighbours_set:
            print(neighbour_id, file=outfile)
    return neighbours_set
        
#neigh_ids = get_bait_neighbourhood("5088749308272be178d0")


def get_p90s(p30_id):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd = "SELECT p90 FROM clusters WHERE p30 = '%s'" % (p30_id)
    cursor.execute(cmd)
    p90s = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(list(set(p90s)))

def get_p100s(p90_id):
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd = "SELECT p100 FROM clusters WHERE p90 = '%s'" % (p90_id)
    cursor.execute(cmd)
    p100s = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(list(set(p100s)))


def calc_icity(pid, neighb_set):
    p30=get_permissive_rep(pid)
    if p30 == None:
        return None # to do - fix this issue by making comprehensive permissive clusters
    p90s = get_p90s(p30)
    #print(p90s)
    hits = 0
    #print("start calculating icity")
    for p90 in p90s:
        p100s = get_p100s(p90)
        #print(p100s)
        for p100 in p100s:
            if p100 in neighb_set:
                hits += 1
                break
    #print("_icity for " + pid + " calculated")
    return [hits / len(p90s), hits, len(p90s)]


# In[30]:


# perm_rep_ex = get_permissive_rep("e0f58eed15ffda8a926c")
# related_baits_ex = get_related_baits(perm_rep_ex)
# bait_neighbourhood_ex = get_full_bait_neighbourhood(related_baits_ex)


# In[31]:


def calc_icity_pool(pid_set_list):
    pool = Pool(cpu_count())
    results = pool.starmap(calc_icity, iterable = pid_set_list)
    pool.close()
    pool.join()
    return results


# In[ ]:





# In[32]:


# def final_icity_output(target_p100ids, bait_p100ids):
#     bait_neighbourhood_ex = set()
#     bait_p30s = []
#     for bait_p100id in bait_p100ids:
#         perm_rep_ex = get_permissive_rep(bait_p100id)
#         bait_p30s.append(perm_rep_ex)
#         related_baits_ex = get_related_baits(perm_rep_ex)
#         bait_neighbourhood_ex.update(get_full_bait_neighbourhood(related_baits_ex))
#     #print(sys.getsizeof(bait_neighbourhood_ex)) # to do - delete later
#     icity_arglist = [(target_p100id, bait_neighbourhood_ex) for target_p100id in target_p100ids]
#     print(target_p100ids)
#     icity_list = calc_icity_pool(icity_arglist)
#     print(len(bait_neighbourhood_ex))
#     with open("rando_prot_icity_output.csv", "w") as outfile:
#         print("target_30id, icity, numer, denom", file=outfile)
#         for i in range(len(icity_arglist)):
#             target_p100id = icity_arglist[i][0]
#             target_p30id = get_permissive_rep(target_p100id)
#             bait_p30ids = str(bait_p30s)
#             #print(icity_list)
#             try:
#                 icity, numer, denom = str(icity_list[i][0]), str(icity_list[i][1]), str(icity_list[i][2])
#                 print(",".join([target_p30id, icity, numer, denom]), file=outfile)
#             except:
#                 pass


# ### tnpb data exploration

# In[ ]:





# In[33]:


tnpBs_list = []
infile_tnpBs = "tnpBs_in_testdb.p100.1e4.txt"
with open(infile_tnpBs, "r") as infile:
    lines = infile.readlines()
    for line in lines:
        p100 = line.split('\n')[0]
        tnpBs_list.append(p100)
len(tnpBs_list)


# In[85]:


def get_permissive_rep_pool(bait_p100s):
    pool = Pool(cpu_count())
    results = pool.map(get_permissive_rep, bait_p100s)
    pool.close()
    pool.join()
    return results


# In[ ]:





# In[86]:


# perm_reps = get_permissive_rep_pool(tnpBs_list[:50])


# In[87]:


def get_related_baits_pool(perm_rep_list):
    pool = Pool(cpu_count())
    results = pool.map(get_related_baits, perm_rep_list)
    pool.close()
    pool.join()
    rv_list = []
    for result in results:
        rv_list += result
    return rv_list


# In[88]:


# related_baits = get_related_baits_pool(perm_reps)
# len(related_baits)


# In[89]:


def get_bait_neighb_pool(related_baits):
    pool = Pool(cpu_count())
    results = pool.map(get_bait_neighbourhood, related_baits)
    pool.close()
    pool.join()
    merged_result_dict = {}
    for result_dict in results:
        merged_result_dict.update(result_dict)
    #print(merged_result_dict)
    return merged_result_dict


# In[90]:


def create_tnpB_neighbourhood(bait_p100ids):
    print("get permissive reps")
    perm_reps = get_permissive_rep_pool(bait_p100ids)
    print("get related baits")
    related_baits = get_related_baits_pool(perm_reps)
    print("get bait neighbs")
    bait_neighbs = get_bait_neighb_pool(related_baits)
    print("writing bait neighbs")
    with open('bait_neighbourhood.tsv', 'w') as outfile:
        for key in bait_neighbs.keys():
            outfile.write("%s\t%s\n"%(key,bait_neighbs[key]))
    # with open("bait_neighbourhood.txt", "w") as outfile:
    #     for neighbour_id in bait_neighbs:
    #         print(neighbour_id, file=outfile)
    return bait_neighbs


# In[91]:


# tic = time.time()
# create_tnpB_neighbourhood(tnpBs_list)
# toc=time.time()


# In[6]:


# required to calculate tnpb icity
# to do - make bait_neighb_ids a globally accessible variable later
# bait_neighb_ids_dict = {}
# with open("bait_neighbourhood.tsv", "r") as infile:
#     lines = infile.readlines()
#     for line in lines:
#         line_list = line.strip("\n").split("\t")
#         neighb_id = line_list[0]
#         bait_ids = ast.literal_eval(line_list[1])
#         bait_neighb_ids_dict[neighb_id] = bait_ids
# bait_neighb_ids = set()
# for bait_neighb_id in bait_neighb_ids_dict.keys():
#     bait_neighb_ids.add(bait_neighb_id)


# In[ ]:


#len(bait_neighb_ids_dict), len(tnpBs_list)


# In[ ]:


# tic = time.time()
# get_baits_closetotarget("133307e8c61ec09a60")
# toc=time.time()


# In[272]:


# get_baits_closetotarget("133307e8c61ec09a60")


# In[93]:


def calc_tnpb_icity(pid, bait_neighb_ids):
    p30 = get_permissive_rep(pid)
    if p30 == None:
        return None # to do - fix this issue by making comprehensive permissive clusters
    p90s = get_p90s(p30)
    hits = 0
    #print("start calculating icity")
    for p90 in p90s:
        p100s = get_p100s(p90)
        for p100 in p100s:
            if p100 in bait_neighb_ids:
                hits += 1
                break
    #print("_icity for " + pid + " calculated")
    return [hits / len(p90s), hits, len(p90s)]
def calc_tnpb_icity_pool(pid_list):
    pool = Pool(cpu_count())
    results = pool.starmap(calc_tnpb_icity, iterable = pid_list)
    pool.close()
    pool.join()
    return results


# In[94]:


def get_baits_closetotarget(p30targetid, bait_neighb_ids_dict):      
    p100s = set()
    p90s = get_p90s(p30targetid)
    for p90 in p90s:
        p100s.update(get_p100s(p90))
    nearby_baits = set()
    for p100 in p100s:
        try:
            baitids = bait_neighb_ids_dict[p100]
            for baitid in baitids:
                nearby_baits.add(baitid)
        except:
            pass
        #baitids = bait_neighb_ids_dict[p100]
    return list(nearby_baits)


# In[95]:


# tnpB_path = "tnpB_icity_output.csv"
# tnpB_df = pd.read_csv(tnpB_path).rename(columns = 
#                                         {" icity": "icity", 
#                                          " numer": "numer", 
#                                          " denom": "denom"}).drop_duplicates()
# #tnpB_df = tnpB_df[tnpB_df["icity"] > .7].sort_values(["icity","numer"], ascending = False)
# target_p30s = list(tnpB_df['target_p100id'])
# len(target_p30s)


# In[106]:


inactive_tnpBs_list = []
with open("../tnpBs/_all_inactive_tnpBs.3.faa", "r") as infile:
    lines=infile.readlines()
    for line in lines:
        if line[0] == '>':
            tnpBid = line.strip('>').strip('\n')
            inactive_tnpBs_list.append(tnpBid)
            
print(len(inactive_tnpBs_list))

with open("../tnpBs/_all_inactive_tnpBs.2.faa", "r") as infile:
    lines=infile.readlines()
    for line in lines:
        if line[0] == '>':
            tnpBid = line.strip('>').strip('\n')
            inactive_tnpBs_list.append(tnpBid)

print(len(inactive_tnpBs_list))

with open("../tnpBs/_all_inactive_tnpBs.1.faa", "r") as infile:
    lines=infile.readlines()
    for line in lines:
        if line[0] == '>':
            tnpBid = line.strip('>').strip('\n')
            inactive_tnpBs_list.append(tnpBid)

print(len(inactive_tnpBs_list))


# In[61]:


test_tnpBs = ['32e56140588692cbfb', '191f1be0b68ca79ff6', '3808fae964963f6310', '318e083e7fb2e55736', 'e88066d62f86ccf767', '2fbae772397b73d65f', 'f3b670d8b2cf9714c7', '6116321d9a36b8933c', '903bc422568dece6a8', '7798f90257f54ece7e', '9acd78fbd6edc88358', '393474d6188f0450c9', 'b18ae9dfef50b1bd8d', '9cd05a18735a3bc637', '8c5a7b9b13d8e56e3a', '6222751607b4b15702', '0145eb10f554ce522d', '81c9c50070bbe1c5cf', '4c43c2cb810631fe3e', '41ef4ef75bd15e3762', '370a9e0aefebf28663', 'bcbe3060b45beaba5f', 'b36f94d2b41b00bf14', '4e5fba4390bd1d67c6', '2b755952ce05675522', 'bebc59761f02536509', 'a2becbc10194ff1ce0']


# In[63]:


len(test_tnpBs), sum([tnpB in test_tnpBs for tnpB in inactive_tnpBs_list])


# In[64]:


sig70_near_dTnpB_3_list = ['aa6af7e9289c3558d3', 'c1050b21cc75640d51', '23422406293d40c201', '6bf6c4c7da68779d7a']


# In[25]:


def tnpB_icity_output(tnpBs_list):
    # maps neighbor genes to nearby baits and bait relatives
    print("getting all genes near baits that could contribute to icity")
    bait_neighb_ids_dict = create_tnpB_neighbourhood(tnpBs_list)
    bait_neighb_ids = set()
    for bait_neighb_id in bait_neighb_ids_dict.keys():
        bait_neighb_ids.add(bait_neighb_id)
    print("mapping all target genes to nearby baits")
    
    ## considers a more restricted pool of tnpB neighbors as target genes
    #target_ids_dict = get_bait_neighb_pool(tnpBs_list)
    target_ids_dict = bait_neighb_ids_dict
    print("target to bait dict created: " + str(len(target_ids_dict)))
    target_p100ids = []
    print(str(len(target_ids_dict.keys())) + " keys in dict")
    for target_id in target_ids_dict.keys():
        target_p100ids.append(target_id)

    icity_arglist = [(target_p100id, bait_neighb_ids) for target_p100id in target_p100ids] 
    print("calculating icity for " + str(len(icity_arglist)) + " target genes")
    icity_list = calc_tnpb_icity_pool(icity_arglist)
    print("icity for " + str(len(icity_list)) + " target genes calculated")
    count_err, getbaits_err = 0, 0
    print("writing -icity to outfile")
    with open("dTnpB_icity_output.tsv", "w") as outfile:
    #with open("sig70_icity_output.tmp.csv", "w") as outfile:
        print("target_p30id\tbaitp100s\ticity\tnumer\tdenom", file=outfile)
        for i in range(len(icity_arglist)):
            target_p100id = icity_arglist[i][0]
            target_p30id = get_permissive_rep(target_p100id)
            try:
                baitp100s = str(get_baits_closetotarget(target_p30id, bait_neighb_ids_dict))
            except:
                getbaits_err += 1
                continue
            try:
                icity, numer, denom = str(icity_list[i][0]), str(icity_list[i][1]), str(icity_list[i][2])
                print("\t".join([target_p30id, baitp100s, icity, numer, denom]), file=outfile)
            except Exception as e:
                count_err += 1
                # print(icity_list[i])
                # print(e)
    print("num getbaits_err is " + str(getbaits_err))
    print("num exceptions is " + str(count_err))
tnpB_icity_output(inactive_tnpBs_list)
#tnpB_icity_output(sig70_near_dTnpB_3_list)


# In[75]:


import numpy as np
import matplotlib.pyplot as plt
plt.hist(sig70_icity_df["baitp100s"].apply(lambda x: x.count(',') + 1))


# In[100]:


sig70_icity_out = "sig70_icity_output.tmp.tsv"
sig70_icity_df = pd.read_csv(sig70_icity_out, sep = "\t")
sig70_icity_df = sig70_icity_df[sig70_icity_df["icity"] > .7]
sig70_icity_df = sig70_icity_df[sig70_icity_df["denom"] > 1]


# In[103]:


len(sig70_icity_df["target_p30id"])


# In[101]:


len(set(sig70_icity_df["target_p30id"]).intersection(tnpBs_list))


# In[102]:


len(set(sig70_icity_df["target_p30id"]).intersection(inactive_tnpBs_list))


# ### get all cas1s and cas2s

# In[ ]:


# path_drep = "../drep_genomes/OUTPUT/rep_genomes/"
# drep_paths = []
# # for directory in os.listdir(path_drep):
# #     drep_samples[directory] = path_drep + "/" + directory
# for directory in os.listdir(path_drep):
#     path_2 = path_drep + directory
#     if os.path.isdir(path_2):
#         for directory_2 in os.listdir(path_2):
#             path_3 = path_2 + "/" + directory_2
#             if os.path.isdir(path_3):
#                 for directory_3 in os.listdir(path_3):
#                     path_4 = path_3 + "/" + directory_3 + "/"
#                     if os.path.isdir(path_4):
#                         for directory_samp in os.listdir(path_4):
#                             samp_dir = path_4 + directory_samp
#                             if os.path.isdir(samp_dir):
#                                 drep_paths.append([samp_dir, directory_samp])


# In[ ]:


def get_castyperfiles(drep_path):
    bucket_path = "gs://durrant/crispestdb/" + drep_path[0][35:] + "/" + drep_path[1] + ".crisprcastyper.domains.tsv.gz"
    command_download = "gsutil cp " + bucket_path + " " + drep_path[0]
    os.system(command_download)
    command_gunzip = "gunzip " + drep_path[0] + "/" + drep_path[1] + ".crisprcastyper.domains.tsv.gz"
    os.system(command_gunzip)
    
def get_castyperfiles_pool(drep_paths):
    pool = Pool(cpu_count())
    results = pool.map(get_castyperfiles, drep_paths)
    pool.close()
    pool.join()

#get_castyperfiles_pool(drep_paths)


# In[ ]:


# drep_cct_paths = []
# for drep_path in drep_paths:
#     cct_path = drep_path[0] + "/" + drep_path[1] + ".crisprcastyper.domains.tsv"
#     drep_cct_paths.append(cct_path)
# len(drep_cct_paths)


# In[ ]:


# # get all cas1s and cas 2s
# count = 0
# for drep_cct_path in drep_cct_paths:
#     try:
#         with open(drep_cct_path, "r") as infile:
#             next(infile)
#             lines = infile.readlines()
#             for line in lines:
#                 line_sep = line.split('\t')
#                 #print(line_sep)
#                 header = line_sep[0]
#                 cas_annot = line_sep[3]
#                 if "Cas1_" in cas_annot:
#                     print(drep_cct_path + "," + header[:20] + "," + cas_annot, file = open("cas1_drep.csv", "a"))
#                 if "Cas2_" in cas_annot:
#                     print(drep_cct_path + "," + header[:20] + "," + cas_annot, file = open("cas2_drep.csv", "a"))
#     except:
#         with open("missed_cct_paths.txt", "a") as outfile:
#             print(drep_cct_path, file=outfile)


# In[ ]:


# cas1s_set, cas2s_set = set(), set()
# with open("cas1_drep.csv", "r") as infile:
#     lines = infile.readlines()
#     for line in lines:
#         line_sep = line.split(',')
#         cas1_id = line_sep[1]
#         cas1s_set.add(cas1_id)
# with open("cas2_drep.csv", "r") as infile:
#     lines = infile.readlines()
#     for line in lines:
#         line_sep = line.split(',')
#         cas2_id = line_sep[1]
#         cas2s_set.add(cas2_id)


# In[ ]:


# len(cas1s_set), len(cas2s_set)


# In[ ]:


def get_random10kprot():
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd = "SELECT hashid FROM proteins ORDER BY hashid LIMIT 10000"
    cursor.execute(cmd)
    random10kprot = get_list_ids_fromcursor(cursor.fetchall())
    conn.close()
    return(random10kprot)


# In[ ]:


# random_9_proteins = ["00000009ba423c2c0a75", "00000051418dfe799d75", "0000005ef429c9ab9d45", "00000071a6b842424f6a", "000000777af16a869fb0", "000000d70aca3d611da1", "000000e13071a7ecf39f", "0000011dab156def9db5", "00000129e4bc5caa1562"]
# final_icity_output(random_9_proteins, list(cas2s_set))


# In[ ]:


def final_icity_output(target_p100ids, bait_p100ids, outfilename):
    bait_neighbourhood_ex = set()
    bait_p30s = []
    for bait_p100id in bait_p100ids:
        perm_rep_ex = get_permissive_rep(bait_p100id)
        bait_p30s.append(perm_rep_ex)
        related_baits_ex = get_related_baits(perm_rep_ex)
        bait_neighbourhood_ex.update(get_full_bait_neighbourhood(related_baits_ex))
    icity_arglist = [(target_p100id, bait_neighbourhood_ex) for target_p100id in target_p100ids]
    icity_list = calc_icity_pool(icity_arglist)
    with open(outfilename, "w") as outfile:
        print("target_30id, bait_30ids, icity, numer, denom", file=outfile)
        for i in range(len(icity_arglist)):
            target_p100id = icity_arglist[i][0]
            target_p30id = get_permissive_rep(target_p100id)
            bait_p30ids = str(bait_p30s)
            asdf
            try:
                icity, numer, denom = str(icity_list[i][0]), str(icity_list[i][1]), str(icity_list[i][2])
                print(",".join([target_p30id, bait_p30ids, icity, numer, denom]), file=outfile)
            except:
                pass


# In[174]:


tnpaid = ['1398176b7a94c33a68']
tnpbid = []
with open('tnpBs_in_testdb.p100.1e4.txt', 'r') as infile:
    lines = infile.readlines()
    for line in lines:
        line = line.strip('\n')
        tnpbid.append(line)
#cas1id = ['6c41d89d162aad350fac']


# In[ ]:




#final_icity_output(tnpbid, tnpaid, "tnpB_tnpA_icity.csv")


# In[175]:


df_highicitysearch = pd.read_csv('tnpB_tnpA_icity.csv')
df_highicitysearch = df_highicitysearch.sort_values([' icity', ' numer'], ascending = False)
df_highicitysearch[df_highicitysearch[' icity'] < 1]


# In[ ]:


# cas1_icity_list = []
# with open("cas2icity_output.csv", "r") as infile:
#     next(infile)
#     lines = infile.readlines()
#     for line in lines:
#         line_sep = line.split(',')
#         cas1_icity = line_sep[1]
#         cas1_icity_list.append(float(cas1_icity))
# cas2_icity_list = []
# with open("cas1icity_output.csv", "r") as infile:
#     next(infile)
#     lines = infile.readlines()
#     for line in lines:
#         line_sep = line.split(',')
#         cas2_icity = line_sep[1]
#         cas2_icity_list.append(float(cas2_icity))
# rando_prot_icity_list = []
# with open("rando_prot_icity_output.csv", "r") as infile:
#     next(infile)
#     lines = infile.readlines()
#     for line in lines:
#         line_sep = line.split(',')
#         rando_prot_icity = line_sep[1]
#         rando_prot_icity_list.append(float(rando_prot_icity))


# In[12]:


# plt.hist(rando_prot_icity_list, bins = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
# plt.xlabel("icity score of random proteins")
# plt.ylabel("frequency")
# plt.title("cas2-icity scores for 10k random proteins in 80k isolate genomes")


# In[13]:


# plt.hist(cas1_icity_list, bins = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
# plt.xlabel("cas2-icity score")
# plt.ylabel("frequency")
# plt.title("cas2-icity scores for all cas1s in 80k isolate genomes")


# In[14]:


# plt.hist(cas2_icity_list, bins = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
# plt.xlabel("cas1-icity score")
# plt.ylabel("frequency")
# plt.title("cas1-icity scores for all cas2s in 80k isolate genomes")


# ### explore icity outputs

# In[ ]:


#final_icity_output(["b9a340f0084638555260"],["e0f58eed15ffda8a926c"])


# In[ ]:


#icity_arglist = [(bait_id, bait_neighbourhood_ex) for bait_id in bait_neighbourhood_ex]


# In[ ]:


# for i in range(10):
#     print(calc_icity(icity_arglist[i][0], icity_arglist[i][1]))
# icity_arglist[0][0]


# In[ ]:


#calc_icity_pool(icity_arglist[:1000])


# ### calculate unrefined p100, p90, p30 -icity

# In[ ]:


def p100icity():
    pass


# In[ ]:




