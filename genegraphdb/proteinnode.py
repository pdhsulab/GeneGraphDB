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
import time
from collections import deque
import csv
from csv import reader

def connect_proteins(coords_csv, max_distance, gene_neighbs = True):
    print("Loading protein2protein edges...")
    tic = time.time()
    outfile = open("protein2protein.tmp.csv", "w")
    print("recid,phash,qhash", file=outfile)
    os.system("cat gene_coords.tmp.csv | sed -e '1s/phash/hash/' | cut -d',' -f 1-5 > gene_coords_m.tmp.csv")
    os.system("cat crispr_coords.tmp.csv | awk 'FNR > 1' > tmp.crispr_coords.csv")
    os.system("cat gene_coords_m.tmp.csv tmp.crispr_coords.csv | sort --field-separator=',' -k1,1 -k4,4n "
              "> merged_sorted_coords.tmp.csv")
    os.system("rm gene_coords_m.tmp.csv tmp.crispr_coords.csv")
    merge_sorted_coords_csv = "merged_sorted_coords.tmp.csv"
    create_protein_pair_csv(merge_sorted_coords_csv, max_distance, outfile, gene_neighbs)
    outfile.close()
    load_csv()
    toc = time.time()
    # remove("protein2protein.tmp.csv")
    # remove(merge_sorted_coords_csv)
    print("Loading protein2protein edges took %f seconds" % (toc - tic))

# saves a temp .csv with all the protein pairs
def create_protein_pair_csv(coords_csv, max_distance, outfile, gene_neighbs):
    with open(coords_csv, 'r') as f:
        infile = reader(f)
        header = next(infile)
        gene_neigh_queue, base_neigh_queue = deque(), deque()
        row1 = next(infile)
        recid, old_phash, old_chash, old_start_coord = row1[0], row1[1], row1[2], row1[3]
        if gene_neighbs:
            gene_neigh_queue.appendleft(old_phash)
            if header != None:
                for line in infile:
                    recid, cur_phash, cur_chash, cur_start_coord = line[0], line[1], line[2], line[3]
                    for old_phash in gene_neigh_queue:
                        print(recid + "," + cur_phash + "," + old_phash, file=outfile)
                    gene_neigh_queue = update_gene_neigh_queue(gene_neigh_queue, cur_phash, old_chash,
                                                      cur_chash, max_distance, cur_start_coord,
                                                      old_start_coord)
                    old_start_coord = cur_start_coord
                    old_chash = cur_chash
        else:
            print("yay! bases!")
            base_neigh_queue.appendleft({"phash": old_phash, "start_coord": old_start_coord})
            if header != None:
                for line in infile:
                    print(line)
                    recid, cur_phash, cur_chash, cur_start_coord = line[0], line[1], line[2], line[3]
                    for _dict in base_neigh_queue:
                        old_phash = _dict["phash"]
                        print(recid + "," + cur_phash + "," + old_phash, file=outfile)
                    base_neigh_queue = update_base_neigh_queue(base_neigh_queue, cur_phash, old_chash,
                                                           cur_chash, max_distance, cur_start_coord,
                                                           old_start_coord)
                    old_start_coord = cur_start_coord
                    old_chash = cur_chash


def update_gene_neigh_queue(queue, cur_phash, old_chash, cur_chash, max_distance, new_coord,
                        old_coord):
    if newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
        if len(queue) >= max_distance:
            queue.pop()
        queue.appendleft(cur_phash)
    else:
        sort_coords_csv()
    return queue

def update_base_neigh_queue(queue, cur_phash, old_chash, cur_chash, max_distance, new_coord,
                        old_coord):
    if newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
        print([old_chash, cur_chash, new_coord, old_coord])
        while len(queue) > 0 and int(new_coord) - int(queue[-1]["start_coord"]) > max_distance:
            queue.pop()
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    else:
        # to do - this is not an exception - handle cases where moving b/w contigs
        print([old_chash, cur_chash, new_coord, old_coord])
        raise Exception("Input gff is not sorted - system bug")
    return queue

# still necessary - checks whether pointer is moving to a new contig
def newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
    try:
        #to do- test if this works
        return int(new_coord) > int(old_coord) or old_chash != cur_chash
    except TypeError:
        print("Type error with gene coordinates")

#csv will have two columns - one for donor protein's phash, other for recipient
def load_csv():
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd_protein2protein_edges = """
              USING PERIODIC COMMIT
              LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
              MATCH (p:Protein), (q:Protein)
              WHERE p.hashid = row.phash AND q.hashid = row.qhash
              MERGE (p)-[f:protein2protein]->(q)
              """.format(
        csv=abspath('protein2protein.tmp.csv')
    )
    cmd_crispr2protein_edges = """
              USING PERIODIC COMMIT
              LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
              MATCH (p:Protein), (c:CRISPR)
              WHERE (p.hashid = row.phash AND c.hashid = row.qhash)
              OR (p.hashid = row.qhash AND c.hashid = row.phash)
              MERGE (p)-[f:protein2crispr]->(c)
              """.format(
        csv=abspath('protein2protein.tmp.csv')
    )
    print(cmd_protein2protein_edges)
    conn.query(cmd_protein2protein_edges, db=DBNAME)
    print(cmd_crispr2protein_edges)
    conn.query(cmd_crispr2protein_edges, db=DBNAME)
    conn.close()

def load_protein_coords():
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = """
          USING PERIODIC COMMIT
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
          MATCH (p:Protein), (c:Contig) 
          WHERE p.hashid = row.phash AND c.hashid = row.chash 
          MERGE (p)-[r:GeneCoord {{start: row.start, end: row.end, orient: row.orient}}]->(c)
          """.format(
        csv=abspath('gene_coords.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()

