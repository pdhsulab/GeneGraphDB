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
from csv import reader


def connect_proteins(coords_csv, max_distance, gene_neighbs = True):
    print("Loading protein2protein edges...")
    tic = time.time()
    outfile = open("protein2protein.tmp.csv", "w")
    print("phash,qhash", file=outfile)
    if gene_neighbs:
        gene_distance_csv(coords_csv, max_distance, outfile)
    else:
        base_distance_csv(coords_csv, max_distance, outfile)
    outfile.close()
    load_csv()
    toc = time.time()
    # remove("protein2protein.tmp.csv")
    print("Loading protein2protein edges took %f seconds" % (toc - tic))

# saves a temp .csv with all the protein pairs
def gene_distance_csv(coords_csv, max_distance, outfile):
    with open(coords_csv, 'r') as f:
        infile = reader(f)
        header = next(infile)

        neigh_queue = deque()
        row1 = next(infile)
        old_phash = row1[0]
        old_chash = row1[1]
        old_start_coord = row1[2]  # only initialize this once, even with many contigs per sample
        print("old_chash is " + str(old_chash))
        print("old_phash is " + str(old_phash))
        neigh_queue.appendleft(old_phash)

        if header != None:
            for line in infile:
                # FIGURE THIS OUT
                cur_phash, cur_chash, cur_start_coord = line[0], line[1], line[2]
                for old_phash in neigh_queue:
                    print(cur_phash + "," + old_phash, file=outfile)
                neigh_queue = update_neighb_queue(neigh_queue, cur_phash, old_chash, cur_chash,
                                                  max_distance, cur_start_coord, old_start_coord)

                old_start_coord = cur_start_coord
                # old_phash = cur_phash
                old_chash = cur_chash

def base_distance_csv(coords_csv, max_distance):
    #TO DO
    pass

def update_neighb_queue(queue, cur_phash, old_chash, cur_chash, max_distance, new_coord, old_coord):
    if newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
        if len(queue) >= max_distance:
            queue.pop()
        queue.appendleft(cur_phash)
    else:
        sort_coords_csv()
    return queue

def newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
    try:
        #to do- test if this works
        return int(new_coord) > int(old_coord) or old_chash != cur_chash
    except:
        print(str(new_coord) + "," + str(old_coord))

def sort_coords_csv():
    # TO DO: sort gene_coords.tmp.csv
    print("sort stuff here")
    pass

#csv will have two columns - one for donor protein's phash, other for recipient
def load_csv():
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd_load_edges = """
              USING PERIODIC COMMIT
              LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
              MATCH (p:Protein), (q:Protein)
              WHERE p.hashid = row.phash AND q.hashid = row.qhash
              MERGE (p)-[f:protein2protein]->(q)
              """.format(
        csv=abspath('protein2protein.tmp.csv')
    )
    print(cmd_load_edges)
    conn.query(cmd_load_edges, db=DBNAME)
    conn.close()
