#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os
from Bio import SeqIO
import pandas as pd


# In[4]:


infile_path = "out/hmmer.tab"
outfile_path = "OUTPUT/596adb7fa02d385b18.examples.gff"


# In[12]:


def parse_cct_hmm(infile_path, outfile_path):
    with open(infile_path, "r") as infile, open(outfile_path, "a") as out_gff:
        next(infile)
        for line in infile:
            line_elem = line.split('\t')
            Hmm, ORF, tlen, qlen, Eval, score, start, end, Acc, Pos, Cov_seq, Cov_hmm, strand = line_elem
            if int(strand) == 1:
                strand = '+'
            elif int(strand) == -1:
                strand = '-'
            seqname, source, feature, score, frame = ORF, "CCT", "CCT_domain", 1, 0
            seqname = '_'.join(seqname.split('_')[:-1])
            attrib = ';'.join(['ID=' + Hmm, 
                              'full_seq_evalue=' + str(Eval), 
                              'full_seq_score=' + str(score),
                              'qlen=' + str(qlen),
                              'tlen=' + str(tlen),
                              'qaln_perc=' + str(float(Cov_hmm)*100),
                              'taln_perc=' + str(float(Cov_seq)*100),
                             ])
            out = [seqname, source, feature, start, end, score, strand, frame, attrib]
            if float(Eval) < 1e-4:
                print(*out, sep = '\t', file=out_gff)


# In[11]:


#parse_hmmsearch('hmmer/out_clean/7205dff302ff900300.csv', 'tmp.gff')


# In[9]:


def parse_hmmsearch(infile_path, outfile_path):
    with open(infile_path, "r") as infile, open(outfile_path, "a") as out_gff:
        next(infile)
        for line in infile:
            line_elem = line.split(',')
            query_id,query_accession,query_len,hit_id,hit_score,hit_evalue,target_length,qstart,qend,tstart,tend = line_elem
            seqname, dnastart, dnaend, strand = hit_id.split(';')
            if strand == "+":
                start = int(dnastart) + 3*int(tstart)
                end = int(dnastart) + 3*int(tend)
            elif strand == "-":
                start = int(dnaend) - 3*int(tend)
                end = int(dnaend) - 3*int(tstart)
            source, feature, score, frame = "hmmsearch", "pfam_domain", 1, 0
            attrib = ';'.join(['ID=' + query_id, 
                              'domain_evalue=' + str(hit_evalue), 
                              'domain_score=' + str(hit_score),
                              'qlen=' + str(query_len),
                              'tlen=' + str(target_length),
                               'qstart=' + str(qstart),
                              'qend=' + str(qend),
                              'tstart=' + str(tstart),
                               'tend=' + str(tend),
                             ])
            out = [seqname, source, feature, start, end, score, strand, frame, attrib]
            if float(hit_evalue) < 1e-4:
                print(*out, sep = '\t', file=out_gff)


# In[6]:


def parse_minced(infile_path, outfile_path):
    os.system("cat {} >> {}".format(infile_path, outfile_path))


# In[40]:


def put_gff_together(targetid):
    cct_out_path = 'cct/OUTPUT/' + targetid + '/hmmer.tab'
    minced_out_path = 'minced/' + targetid + '.minced.gff'
    gff_path = 'OUTPUT/' + targetid + '.examples.gff'
    sorted_path = 'OUTPUT/' + targetid + '.sorted.gff'
    fna_path = 'INPUT/fna/' + targetid + '.examples.fna'
    fna_sorted_path = 'INPUT/fna/' + targetid + '.sorted.fna'
    gff_path_final = 'OUTPUT/' + targetid + '.final.gff'
    hmm_outclean_path = 'hmmer/out_clean/' + targetid + '.csv'
    gff_hmmsearch_path = 'hmmer/out_clean/' + targetid + '.gff'
    
    parse_cct_hmm(cct_out_path, gff_path)
    parse_minced(minced_out_path, gff_path)
    parse_hmmsearch(hmm_outclean_path, gff_path)
    
    contigs = {}
    with open(fna_path) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            header = rec.id
            sequence = rec.seq
            contigs[header] = sequence
    headers = list(contigs.keys())
    headers.sort()
    os.system("sortBed -i " + gff_path + " > " + sorted_path)
    with open(fna_sorted_path, 'w') as outfile:
        for header in headers:
            print('>' + header, file=outfile)
            print(contigs[header], file=outfile)
    os.system("cat " + sorted_path + " " + fna_sorted_path + " > " + gff_path_final)
    


# In[33]:


def gff_output_raw_csv(targetid):
    gff_path_final = 'OUTPUT/' + targetid + '.final.gff'
    #to do - finish this summaryt stat csv
    _csv_raw_out = 'OUTPUT/_csv_raw_out.csv' 
    with open(gff_path_final, "r") as infile, open(_csv_raw_out, "w") as outfile:
        lines = infile.readlines()
        for line in lines:
            if line[0] == '>':
                break
            line = line.split('\t')
            ID = line[-1]
            print(ID)
            asdf
            print(line, file = outfile)


# In[ ]:




