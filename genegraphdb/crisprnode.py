from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import hashlib
import time
from collections import deque, ChainMap
from csv import reader
from BCBio import GFF

def load_CRISPRs():
    # create fasta from minced.gff
    print("Loading CRISPRs...")
    tic = time.time()
    outfile_crispr = open("CRISPRs.tmp.csv", "w")
    #outfile_minced = open("minced_hashed.tmp.gff", "w")
    print("hashid,repeat,repeat_len,array_len,num_spacers", file=outfile_crispr)
    done = set()
    with open("temp.minced.gff", "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            line = line.strip().split('\t')
            if len(line) != 9:
                print(len(line))
                continue
            array_repeat = re.findall(r"=(.*)", line[8].split(";")[3])[0]
            repeat_len = len(array_repeat)
            array_len = abs(int(line[4]) - int(line[3])) + 1
            num_spacers = line[5]
            unique_str = array_repeat + "," + str(repeat_len) + "," + str(array_len) + \
                         "," + str(num_spacers)
            crhash = hashlib.sha256(unique_str.encode()).hexdigest()
            name_truncate = crhash[:20]
            unique_str = name_truncate + "," + unique_str
            if unique_str in done:
                continue
            done.add(unique_str)  # is this necessary?
            print(unique_str, file=outfile_crispr)
    outfile_crispr.close()
    del done
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:CRISPR {{hashid: row.hashid, repeat: row.repeat, " \
          "repeat_len: row.repeat_len, array_len: row.array_len, " \
          "num_spacers: row.num_spacers}})".format(csv=abspath('CRISPRs.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()
    toc = time.time()
    #remove('CRISPRs.tmp.csv')
    print("Loading CRISPRs took %f seconds" % (toc-tic))

def merge_gff(sample_id):
    # parse through minced.gff
    # cat 8156401/8156401.minced.gff | grep ID=CRISPR > temp.minced.gff
    # cat temp.prodigal.gff temp.minced.gff > temp.merged.gff
    # sortBed -i temp.merged.gff > temp.merged.sorted.gff
    print("start merging gffs")
    protein_path = str(sample_id) + ".prodigal.gff"
    minced_gff_path = str(sample_id) + ".minced.gff"
    os.system("gunzip -d -c " + protein_path + ".gz > " + protein_path)
    os.system("gunzip -d -c " + minced_gff_path + ".gz > " + minced_gff_path)
    os.system("cat " + minced_gff_path + " | grep ID=CRISPR > temp.minced.gff")
    os.system("cat " + protein_path + " temp.minced.gff > temp.merged.gff")
    os.system("sortBed -i temp.merged.gff > temp.merged.sorted.gff")
    os.system("rm temp.merged.gff " + protein_path + " " + minced_gff_path)
    # os.system("rm temp.merged.sorted.gff")
    print("finished merging gffs")
    return("temp.merged.sorted.gff")