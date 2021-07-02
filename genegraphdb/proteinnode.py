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
from collections import deque, ChainMap
from csv import reader

def connect_proteins(coords_csv, max_distance, gene_neighbs = True):
    print("Loading protein2protein edges...")
    tic = time.time()
    outfile = open("protein2protein.tmp.csv", "w")
    print("phash,qhash", file=outfile)
    create_protein_pair_csv(coords_csv, max_distance, outfile, gene_neighbs)
    outfile.close()
    load_csv()
    toc = time.time()
    # remove("protein2protein.tmp.csv")
    print("Loading protein2protein edges took %f seconds" % (toc - tic))

# saves a temp .csv with all the protein pairs
def create_protein_pair_csv(coords_csv, max_distance, outfile, gene_neighbs):
    with open(coords_csv, 'r') as f:
        infile = reader(f)
        header = next(infile)
        gene_neigh_queue, base_neigh_queue = deque(), deque()
        row1 = next(infile)
        old_phash, old_chash, old_start_coord = row1[0], row1[1], row1[2]
        if gene_neighbs:
            gene_neigh_queue.appendleft(old_phash)
            if header != None:
                for line in infile:
                    cur_phash, cur_chash, cur_start_coord = line[0], line[1], line[2]
                    for old_phash in gene_neigh_queue:
                        print(cur_phash + "," + old_phash, file=outfile)
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
                    cur_phash, cur_chash, cur_start_coord = line[0], line[1], line[2]
                    for _dict in base_neigh_queue:
                        old_phash = _dict["phash"]
                        print(cur_phash + "," + old_phash, file=outfile)
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
        try:
            while int(new_coord) - int(queue[-1]["start_coord"]) > max_distance:
                queue.pop()
        except TypeError:
            print("Type error with gene coordinates")
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    else:
        sort_coords_csv()
    return queue

def newGene_is_sorted(old_chash, cur_chash, new_coord, old_coord):
    try:
        #to do- test if this works
        return int(new_coord) > int(old_coord) or old_chash != cur_chash
    except TypeError:
        print("Type error with gene coordinates")

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
