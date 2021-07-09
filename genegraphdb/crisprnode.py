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

def load_CRISPRs(sample_id):
    # create fasta from minced.gff
    print("Loading CRISPRs...")
    tic = time.time()
    outfile_crispr = open("CRISPRs.tmp.csv", "w")
    print("hashid,repeat_len,array_len,num_spacers", file=outfile_crispr)
    done = set()
    crisprid_to_crhash = dict()
    with open("temp.minced.gff") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            line = line.strip().split('\t')
            if len(line) != 9:
                print(len(line))
                continue
            array_repeat = re.findall(r"=(.*)", line[8].split(";")[3])[0]
            crhash = hashlib.sha256(array_repeat.encode()).hexdigest()
            name_truncate = crhash[:20]
            repeat_len = len(array_repeat)
            array_len = abs(int(line[4]) - int(line[3])) + 1
            num_spacers = line[5]
            unique_str = name_truncate + "," + str(repeat_len) + "," + str(array_len) + \
                         "," + str(num_spacers)
            crisprid = re.findall(r"=(.*)", line[8].split(";")[0])[0]
            crisprid_to_crhash[crisprid] = name_truncate
            if unique_str in done:
                continue
            done.add(unique_str)  # is this necessary?
            print(unique_str, file=outfile_crispr)
    outfile_crispr.close()
    del done
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:CRISPR {{hashid: row.hashid, " \
          "repeat_len: row.repeat_len, array_len: row.array_len, " \
          "num_spacers: row.num_spacers}})".format(csv=abspath('CRISPRs.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()
    toc = time.time()
    os.system('rm CRISPRs.tmp.csv temp.minced.gff')
    print("Loading CRISPRs took %f seconds" % (toc-tic))
    return crisprid_to_crhash

def merge_gff(sample_id):
    # parse through minced.gff
    # cat 8156401/8156401.minced.gff | grep ID=CRISPR > temp.minced.gff
    # cat temp.prodigal.gff temp.minced.gff > temp.merged.gff
    # sortBed -i temp.merged.gff > temp.merged.sorted.gff
    print("start merging gffs")
    protein_path = str(sample_id) + ".prodigal.gff"
    os.system("gunzip -d -c " + protein_path + ".gz > " + protein_path)
    minced_gff_path = str(sample_id) + ".minced.gff"
    os.system("gunzip -d -c " + protein_path + ".gz > " + protein_path)
    os.system("gunzip -d -c " + minced_gff_path + ".gz > " + minced_gff_path)
    os.system("cat " + minced_gff_path + " | grep ID=CRISPR > temp.minced.gff")
    os.system("cat " + protein_path + " temp.minced.gff > temp.merged.gff")
    return_filename = "merged.sorted.tmp.gff"
    os.system("sortBed -i temp.merged.gff > " + return_filename)
    os.system("rm temp.merged.gff " + protein_path)
    # os.system("rm " + return_filename)
    print("finished merging gffs")
    return(return_filename)

def load_crispr_coords():
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = """
          USING PERIODIC COMMIT
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
          MATCH (cr:CRISPR), (c:Contig) 
          WHERE cr.hashid = row.crhash AND c.hashid = row.chash 
          MERGE (cr)-[r:CrisprCoord {{start: row.start, end: row.end}}]->(c)
          """.format(
        csv=abspath('crispr_coords.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()