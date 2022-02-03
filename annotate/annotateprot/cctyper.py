#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import os
from multiprocessing import Pool


# %%


infile_paths = ['INPUT/fna/' + file for file in os.listdir('INPUT/fna') if file != '.ipynb_checkpoints']
p30s = [infile_path.split('/')[-1].split('.')[0] for infile_path in infile_paths]

# %%

def cctyper(infile_paths):
    for infile_path in infile_paths:
        p30 = infile_path.split('/')[-1].split('.')[0]
        outfile_path = 'cct/OUTPUT/' + p30
        cmd = "cctyper " + infile_path + " " + outfile_path + " --prodigal meta"
        os.system(cmd)
cctyper(infile_paths)
