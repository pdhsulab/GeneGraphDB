#!/usr/bin/env python
# coding: utf-8

# In[21]:


import os
from Bio import SeqIO


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
            #print(*out, sep = '\t')
            print(*out, sep = '\t', file=out_gff)


# In[ ]:


def parse_hmmer():
    pass


# In[40]:


def put_gff_together(targetid):
    cct_out_path = 'cct/OUTPUT/' + targetid + '/hmmer.tab'
    gff_path = 'OUTPUT/' + targetid + '.examples.gff'
    sorted_path = 'OUTPUT/' + targetid + '.sorted.gff'
    fna_path = 'INPUT/fna/' + targetid + '.examples.fna'
    fna_sorted_path = 'INPUT/fna/' + targetid + '.sorted.fna'
    gff_path_final = 'OUTPUT/' + targetid + '.final.gff'
    
    parse_cct_hmm(cct_out_path, gff_path)
    
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
    


# In[ ]:




