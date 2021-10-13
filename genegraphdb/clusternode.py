from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from genegraphdb import testing
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

def load_cluster_nodes():
    print("Loading cluster nodes...")
    tic = time.time()
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)

    cmd90 = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:cluster90 {{hashid: row.p90}})".format(
        csv=abspath("clusters_20Krows.csv")
    )
    print(cmd90)
    conn.query(cmd90, db=DBNAME)

    cmd30 = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:cluster30 {{hashid: row.p30}})".format(
        csv=abspath("clusters_20Krows.csv")
    )
    print(cmd30)
    conn.query(cmd30, db=DBNAME)

    cmd_connect_prot_90 = """
        LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
        MATCH (p:Protein), (clust:cluster90) 
        WHERE (p.hashid = row.p100 AND clust.hashid = row.p90)
        MERGE (p)-[r:prot2cluster90]->(clust)
    """.format(
        csv=abspath("clusters_20Krows.csv")
    )
    print(cmd_connect_prot_90)
    conn.query(cmd_connect_prot_90, db=DBNAME)

    cmd_connect_30_90 = """
        LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
        MATCH (c30:cluster30), (c90:cluster90) 
        WHERE (c30.hashid = row.p30 AND c90.hashid = row.p90)
        MERGE (c30)-[r:cluster30tocluster90]->(c90)
    """.format(
        csv=abspath("clusters_20Krows.csv")
    )
    print(cmd_connect_30_90)
    conn.query(cmd_connect_30_90, db=DBNAME)

    conn.close()
    toc = time.time()
    print("Loading cluster nodes took %f seconds" % (toc-tic))