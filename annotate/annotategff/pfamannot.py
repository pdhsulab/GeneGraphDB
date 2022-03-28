#!/usr/bin/env python
# coding: utf-8

# In[8]:


from multiprocessing import Pool, cpu_count
from Bio import SearchIO
import os


# In[41]:


def hmmsearch(tnpB_multifaa_path, out_raw_path, hmmfile_path):
    pfam_command = "hmmsearch --cpu {threads} -E 1 --domtblout {outtbl} {hmmfile} {seqdb}".format(
    threads = cpu_count(), outtbl=out_raw_path, hmmfile=hmmfile_path, seqdb=tnpB_multifaa_path)
    print("running hmmsearch")
    os.system(pfam_command)
def format_hmmsearch_output(pid):
    out_raw_path = 'hmmer/out_raw/' + pid + '.tsv'
    out_clean_path = 'hmmer/out_clean/' + pid + '.csv'
    print("formatting hmmsearch results")
    with open(out_raw_path
              , 'r') as infile, open(out_clean_path, "w") as outfile:
        print("query_id,query_accession,query_len,hit_id,hit_score,hit_evalue,target_length,qstart,qend,tstart,tend", file=outfile)
        for qresult in SearchIO.parse(infile, 'hmmscan3-domtab'):
            query_id = qresult.id  #sequence ID from fasta
            query_accession = qresult.accession
            query_len = qresult.seq_len
            hits = qresult.hits
            num_hits = len(hits)
            if num_hits > 0:
                for hit in hits: 
                    hit_evalue = hit.evalue
                    hit_score = hit.bitscore
                    target_length = hit.seq_len
                    hit_id = hit.id
                    for HSP in hit:
                        out_row = ",".join([query_id, query_accession, str(query_len), hit_id, str(hit_score), str(hit_evalue), str(target_length), str(HSP.hit_start), str(HSP.hit_end), str(HSP.query_start), str(HSP.query_end)])
                        print(out_row, file=outfile)


# In[11]:


tnpB_multifaa_path = 'INPUT/faa/7205dff302ff900300.examples.faa'
out_raw_path = 'hmmer/out_raw/7205dff302ff900300.tsv'
out_clean_path = 'hmmer/out_clean/7205dff302ff900300.csv'
hmmfile_path = '../../../tnpBs/data/Pfam-A.hmm'
#hmmsearch(tnpB_multifaa_path, out_raw_path, hmmfile_path)


# In[ ]:





# In[12]:


def hmmsearch_pool(ids_list):
    paths_list = []
    for targetid in ids_list:
        tnpB_multifaa_path = 'INPUT/faa/' + targetid + '.examples.faa'
        out_raw_path = 'hmmer/out_raw/' + targetid + '.tsv'
        hmmfile_path = '../../../tnpBs/data/Pfam-A.hmm'
        paths_list.append((tnpB_multifaa_path, out_raw_path, hmmfile_path))
    pool = Pool(cpu_count())
    results = pool.starmap(hmmsearch, iterable = paths_list)
    pool.close()
    pool.join()
    return results


# In[ ]:





# In[13]:


ids_list  = ['7205dff302ff900300', '6c3b2ade31d3833745', '23422406293d40c201',
       'c1050b21cc75640d51', 'aa6af7e9289c3558d3', 'd2246f26c16fb9eec9',
       '6bf6c4c7da68779d7a']


# In[ ]:





# In[14]:


#hmmsearch_pool(ids_list)


# In[ ]:




