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

def connect_proteins_crisprs(sample_id, max_distance, gene_neighbs = True):
    sample_id_path = sample_id + "/"
    print("Loading protein2protein edges...")
    tic = time.time()
    outfile = open(sample_id_path + "protein2protein.tmp.csv", "w")
    print("recid,phash,qhash", file=outfile)
    os.system("cat " + sample_id_path + "gene_coords.tmp.csv | sed -e '1s/phash/hash/' | cut -d',' -f 1-5 > "
              + sample_id_path + "gene_coords_m.tmp.csv")
    os.system("cat " + sample_id_path + "crispr_coords.tmp.csv | awk 'FNR > 1' > " + sample_id_path + "tmp.crispr_coords.csv")
    os.system("cat " + sample_id_path + "gene_coords_m.tmp.csv " + sample_id_path + "tmp.crispr_coords.csv > "
              + sample_id_path + "merged_coords.tmp.csv")
    os.system("head -n1 " + sample_id_path + "merged_coords.tmp.csv > " + sample_id_path +
              "merged_sorted_coords.tmp.csv && tail -n+2 " + sample_id_path + "merged_coords.tmp.csv | sort "
              "--field-separator=',' -k1,1 -k4,4n >> " + sample_id_path + "merged_sorted_coords.tmp.csv")
    os.system("rm " + sample_id_path + "gene_coords_m.tmp.csv " + sample_id_path + "tmp.crispr_coords.csv "
              + sample_id_path + "merged_coords.tmp.csv")
    merge_sorted_coords_csv = sample_id_path + "merged_sorted_coords.tmp.csv"
    create_protein_pair_csv(merge_sorted_coords_csv, max_distance, outfile, gene_neighbs)
    outfile.close()
    load_csv(sample_id_path + "protein2protein.tmp.csv")
    toc = time.time()
    # remove("protein2protein.tmp.csv")
    # remove(merge_sorted_coords_csv)
    print("Loading protein2protein edges took %f seconds" % (toc - tic))
    return toc-tic

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
                                                      cur_chash, 1, cur_start_coord,
                                                      old_start_coord)
                    old_start_coord = cur_start_coord
                    old_chash = cur_chash
        else:
            print("yay! bases!")
            base_neigh_queue.appendleft({"phash": old_phash, "start_coord": old_start_coord})
            if header != None:
                for line in infile:
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
    if newGene_is_same_contig(old_chash, cur_chash, new_coord, old_coord):
        if len(queue) >= max_distance:
            queue.pop()
        queue.appendleft(cur_phash)
    else:
        queue = deque()
        queue.appendleft(cur_phash)
    return queue

def update_base_neigh_queue(queue, cur_phash, old_chash, cur_chash, max_distance, new_coord,
                        old_coord):
    if newGene_is_same_contig(old_chash, cur_chash, new_coord, old_coord):
        while (int(new_coord) - int(queue[-1]["start_coord"])) > max_distance:
            queue.pop()
            if len(queue) == 0:
                break
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    else:
        queue = deque()
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    return queue

# still necessary - checks whether pointer is moving to a new contig
def newGene_is_same_contig(old_chash, cur_chash, new_coord, old_coord):
        return old_chash == cur_chash

#csv will have two columns - one for donor protein's phash, other for recipient
def load_csv(csv_path):
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd_protein2protein_edges = """
              USING PERIODIC COMMIT
              LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
              MATCH (p:Protein), (q:Protein)
              WHERE p.hashid = row.phash AND q.hashid = row.qhash
              MERGE (p)-[f:protein2protein]->(q)
              """.format(
        csv=abspath(csv_path)
    )
    cmd_crispr2protein_edges = """
              USING PERIODIC COMMIT
              LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
              MATCH (p:Protein), (c:CRISPR)
              WHERE (p.hashid = row.phash AND c.hashid = row.qhash)
              OR (p.hashid = row.qhash AND c.hashid = row.phash)
              MERGE (p)-[f:protein2crispr]->(c)
              """.format(
        csv=abspath(csv_path)
    )
    print(cmd_protein2protein_edges)
    conn.query(cmd_protein2protein_edges, db=DBNAME)
    print(cmd_crispr2protein_edges)
    conn.query(cmd_crispr2protein_edges, db=DBNAME)
    conn.close()

def load_protein_coords(sample_id):
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = """
          USING PERIODIC COMMIT
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
          MATCH (p:Protein), (c:Contig) 
          WHERE p.hashid = row.phash AND c.hashid = row.chash 
          MERGE (p)-[r:GeneCoord {{start: row.start, end: row.end, orient: row.orient}}]->(c)
          """.format(
        csv=abspath(sample_id + '/gene_coords.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()

