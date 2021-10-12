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
    cmd90 = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:90cluster {{hashid: row.p90}})".format(
        csv=abspath("clusters_20Krows.csv")
    )

    print(cmd90)
    conn.query(cmd90, db=DBNAME)
    cmd30 = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:30cluster {{hashid: row.p30}})".format(
        csv=abspath("clusters_20Krows.csv")
    )
    print(cmd30)
    conn.query(cmd30, db=DBNAME)
    conn.close()
    toc = time.time()
    remove(sample_id + '/proteins.tmp.csv')
    print("Loading proteins took %f seconds" % (toc-tic))