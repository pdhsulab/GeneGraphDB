#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.Seq import Seq
import os


# In[2]:


infile_paths = ['INPUT/gff/' + file for file in os.listdir('INPUT/gff') if file != '.ipynb_checkpoints']
infile_paths


# In[7]:


# integrate later
def single_gff_to_fna(infile_path):
    multifasta_path = infile_path.replace("gff", "fna")
    start_reading = False
    outfile_gff_path = infile_path.replace("INPUT/gff", "OUTPUT")
    with open(infile_path, 'r') as infile, open(outfile_gff_path, 'w') as outfile_gff, open(multifasta_path, "w") as outfile:
        for line in infile:
            line = line.strip('\n')
            if line.startswith(">"):
                start_reading = True
            if not start_reading:
                print(line, file=outfile_gff)
            if start_reading:
                print(line, file=outfile)


# In[67]:


def single_gff_to_faa(infile_path):
    faa_path = infile_path.replace("gff", "faa")
    fna_path = infile_path.replace("gff", "fna")
    contigs = fasta_id_dict(fna_path)
    #outfile_gff_path = infile_path.replace("INPUT/gff", "OUTPUT")
    with open(infile_path, 'r') as infile, open(faa_path, "w") as outfile:
        for line in infile:
            line = line.strip('\n').split('\t')
            if 'PRODIGAL' in line:
                header_old = line[0]
                header = ';'.join([line[0], line[3], line[4], line[6]])
                strand = line[6]
                try:
                    dna_seq = contigs[header_old]
                except:
                    pass #some proteins in contigs are annotated without reference sequence
                aa_seq = extract_faa_seq(header, dna_seq)
                print('>' + header, file=outfile)
                print(aa_seq, file=outfile)
                


# In[68]:


def extract_faa_seq(header, seq):
    header_list = header.split(';')
    start = int(header_list[-3]) - 1
    end = int(header_list[-2])
    strand = header_list[-1]
    if strand == '+':
        dna = Seq(seq[start:end])
        protein = dna.translate()
        return protein
    elif strand == '-':
        dna = Seq(seq[start:end]).reverse_complement()
        protein = dna.translate()
        return protein
    
def fasta_id_dict(fna_path):
    contigs = {}
    with open(fna_path) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            header = rec.id
            sequence = rec.seq
            contigs[header] = sequence
    return contigs
        
    


# In[ ]:




