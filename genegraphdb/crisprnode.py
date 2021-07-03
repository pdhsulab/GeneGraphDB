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

def load_CRISPRs(sample_id, crispr):
    # create fasta from minced.gff
    print("Loading CRISPRs...")
    tic = time.time()
    outfile = open("CRISPRs.tmp.csv", "w")

    minced_gff_path = str(sample_id) + ".minced.gff.gz"
    os.system("cat " + minced_gff_path + " | grep ID=CRISPR > temp.minced.gff.gz")

    print("hashid,repeat,repeat_len,array_len,num_spacers", file=outfile)
    done = set()
    os.system("gunzip temp.minced.gff.gz")
    in_handle = open("temp.minced.gff")
    out_handle = open("temp.minced.gff")

    #TO DO - write new hash IDs to out_handle, temp.minced.gff
    for rec in GFF.parse(in_handle):
        for cr_index in range(len(rec.features)):
            array_repeat = rec.features[cr_index].qualifiers["rpt_unit_seq"][0]
            crhash = hashlib.sha256(array_repeat.encode()).hexdigest()
            rec.features[cr_index].qualifiers["ID"] = crhash
            #print(rec.features[cr_index].qualifiers["ID"])
            name_truncate = crhash[:20]

            repeat_len = len(array_repeat)
            array_coords_str = rec.features[cr_index].location
            array_len = abs(array_coords_str.start - array_coords_str.end)
            num_spacers = int(rec.features[cr_index].qualifiers["score"][0])
            unique_str = name_truncate+","+array_repeat+","+str(repeat_len)+","+str(array_len)+\
                         ","+str(num_spacers)
            if unique_str in done:
                continue
            done.add(unique_str) #is this necessary?
            print(unique_str, file=outfile)
    outfile.close()
    del done
    in_handle.close()
    # conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    # cmd = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:Protein {{hashid: row.hashid, length:row.length}})".format(
    #     csv=abspath('proteins.tmp.csv')
    # )
    # print(cmd)
    # conn.query(cmd, db=DBNAME)
    # conn.close()
    # toc = time.time()
    # remove('CRISPRs.tmp.csv')
    #print("Loading CRISPRs took %f seconds" % (toc-tic))

def merge_gff(sample_id):
    # parse through minced.gff
    # cat 8156401/8156401.minced.gff | grep ID=CRISPR > temp.minced.gff
    # cat temp.prodigal.gff temp.minced.gff > temp2.merged.gff
    # sortBed -i temp.merged.gff > temp.merged.sorted.gff
    print("start merging gffs")
    protein_path = str(sample_id) + ".prodigal.gff.gz"
    # minced_gff_path = str(sample_id) + ".minced.gff"
    # os.system("cat " + minced_gff_path + " | grep ID=CRISPR > temp.minced.gff")
    os.system("cat " + protein_path + " temp.minced.gff.gz > temp.merged.gff.gz")
    # unzip?
    os.system("gunzip temp.merged.gff.gz")
    os.system("sortBed -i temp.merged.gff > temp.merged.sorted.gff")
    # os.system("rm temp.minced.gff.gz temp.merged.gff.gz temp.merged.gff")
    # os.system("rm temp.merged.sorted.gff")
    print("finished merging gffs")