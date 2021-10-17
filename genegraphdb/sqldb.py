from genegraphdb import *
from genegraphdb import graphdb
import os
import subprocess
import sqlite3

def hasdb():
    if "genegraph.db" in os.listdir():
        return True
    return False

def createdb():
    print("Creating database: %s" % DBNAME)
    con = sqlite3.connect('genegraph.db')
    cur = con.cursor()
    cur.execute('''CREATE TABLE proteins (hashid text, length real, PRIMARY KEY(hashid, length))''')
    cur.execute('''CREATE TABLE samples (sampleid text, PRIMARY KEY (sampleid))''')
    cur.execute('''CREATE TABLE crisprs (hashid text, PRIMARY KEY (hashid))''')
    cur.execute('''CREATE TABLE contigs (hashid text, length real, PRIMARY KEY (hashid, length))''')
    cur.execute('''CREATE TABLE contig2sample (contighashid text, sampleid text, PRIMARY KEY (contighashid, sampleid))''')
    cur.execute('''CREATE TABLE crisprcoords (crisprhash text, contighash text, start real, end real, PRIMARY KEY (crisprhash, contighash, start))''')
    cur.execute('''CREATE TABLE proteincoords (phash text, contighash text, start real, end real, orientation text, PRIMARY KEY (phash, contighash, start))''')
    cur.execute('''CREATE TABLE prot2prot (p1hash text, p2hash text, PRIMARY KEY (p1hash, p2hash))''')
    cur.execute('''CREATE TABLE prot2crispr (p1hash text, crisprhash text, PRIMARY KEY (p1hash, crisprhash))''')
    cur.execute('''CREATE TABLE prot2protwindow (p1hash text, p2hash text, PRIMARY KEY (p1hash, p2hash))''')
    cur.execute('''CREATE TABLE prot2crisprwindow (p1hash text, crisprhash text, PRIMARY KEY (p1hash, crisprhash))''')
    con.close()
    #os.system("sqlite3 genegraph.db")

    print("SQL database created")

def cleardb():
    os.system("rm genegraph.db")
    pass

# def num_nodes():
#     con = sqlite3.connect('genegraph.db')
#     cur = con.cursor()
#     num_nodes = cur.execute() #conn.query("MATCH (n) RETURN count(*)", db=DBNAME)[0]['count(*)']
#     con.close()
#
#     return num_nodes
#
# def num_rels():
#     con = sqlite3.connect('genegraph.db')
#     cur = con.cursor()
#     num_rels = cur.execute() #conn.query("MATCH ()-[r]->() RETURN count(r) as count", db=DBNAME)[0]['count']
#     con.close()
#
#     return num_rels