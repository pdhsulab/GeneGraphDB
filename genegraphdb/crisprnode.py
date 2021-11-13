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

def merge_gff(sample_id, samples_path):
    # parse through minced.gff
    # cat 8156401/8156401.minced.gff | grep ID=CRISPR > temp.minced.gff
    # cat temp.prodigal.gff temp.minced.gff > temp.merged.gff
    # sortBed -i temp.merged.gff > temp.merged.sorted.gff
    sample_id_path = sample_id + "/"
    print("start merging gffs")
    protein_path = samples_path + sample_id_path + str(sample_id) + ".prodigal.gff"
    os.system("gunzip -d -c " + protein_path + ".gz > " + protein_path)
    minced_gff_path = samples_path + sample_id_path + str(sample_id) + ".minced.gff"
    os.system("gunzip -d -c " + protein_path + ".gz > " + protein_path)
    os.system("gunzip -d -c " + minced_gff_path + ".gz > " + minced_gff_path) # TO DO - change this file name?
    os.system("cat " + minced_gff_path + " | grep ID=CRISPR > " + samples_path + sample_id_path + "temp.minced.gff")
    os.system("cat " + protein_path + " " + samples_path + sample_id_path + "temp.minced.gff > " + samples_path + sample_id_path + "temp.merged.gff")
    return_filename = samples_path + sample_id_path + "merged.sorted.tmp.gff"
    os.system("sortBed -i " + samples_path + sample_id_path + "temp.merged.gff > " + return_filename)
    os.system("rm " + samples_path + sample_id_path + "temp.merged.gff " + protein_path)
    # os.system("rm " + return_filename)
    print("finished merging gffs")
    return(return_filename)

def load_CRISPRs(sample_id, samples_path):
    # create fasta from minced.gff
    print("Loading CRISPRs...")
    tic = time.time()
    outfile_crispr = open(samples_path + sample_id + "/CRISPRs.tmp.csv", "w")
    #print("hashid,repeat_len,array_len,num_spacers", file=outfile_crispr)
    print("hashid", file=outfile_crispr)
    done = set()
    crisprid_to_crhash = dict()
    with open(samples_path + sample_id + "/temp.minced.gff") as infile:
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
            # if unique_str in done:
            #     continue
            # done.add(unique_str)  # is this necessary?
            if name_truncate in done:
                continue
            done.add(name_truncate)
            # print(unique_str, file=outfile_crispr)
            print(name_truncate, file=outfile_crispr)
    outfile_crispr.close()
    del done
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    # to do - make crispr nodes simpler? only include coordinates and pointers to contigs as appropriate
    # cmd = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:CRISPR {{hashid: row.hashid, " \
    #       "repeat_len: row.repeat_len, array_len: row.array_len, " \
    #       "num_spacers: row.num_spacers}})".format(csv=abspath(sample_id + '/CRISPRs.tmp.csv')
    # )
    #conn.execute("PRAGMA journal_mode=WAL")
    cmd = """
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
          MERGE (n:CRISPR {{hashid: row.hashid}})
          """.format(csv=abspath(samples_path + sample_id + '/CRISPRs.tmp.csv'))
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()
    toc = time.time()
    # to do - uncomment this rm command below
    # os.system('rm ' + sample_id + '/CRISPRs.tmp.csv ' + sample_id + '/temp.minced.gff')
    print("Loading CRISPRs took %f seconds" % (toc-tic))
    return crisprid_to_crhash

def load_crispr_coords(sample_id, samples_path):
    sample_id_path = samples_path + sample_id + "/"
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    #conn.execute("PRAGMA journal_mode=WAL")
    cmd = """
          USING PERIODIC COMMIT
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
          MATCH (cr:CRISPR), (c:Contig) 
          WHERE (cr.hashid = row.crhash AND c.hashid = row.chash)
          MERGE (cr)-[r:CrisprCoord {{start: row.start, end: row.end}}]->(c)
          """.format(
        csv=abspath(sample_id_path + 'crispr_coords.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()