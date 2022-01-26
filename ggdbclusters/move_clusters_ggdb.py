#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time


# In[ ]:


def move_clusters_to_genegraphdb():
    conn = sqlite3.connect('genegraph.db')
    cursor = conn.cursor()
    cmd1 = "DROP TABLE stringent"
    cmd2 = "ATTACH DATABASE 'clusters.db' AS db2"
    #cmd3 = "CREATE TABLE permissive (reppid text, pid text, PRIMARY KEY(pid))"
    #cmd4 = "INSERT INTO permissive SELECT a.reppid, a.pid FROM db2.permissive a"
    cmd5 = "CREATE TABLE stringent (reppid text, pid text, PRIMARY KEY(pid))"
    cmd6 = "INSERT INTO stringent SELECT SUBSTR(b.reppid, 1, 18), SUBSTR(b.pid, 1, 18) FROM db2.stringent_ b"
    #cmd7 = "CREATE TABLE clusters (p30 text, p90 text, p100 text, PRIMARY KEY(p100))"

    clusters_csv = open("../clusters/OUTPUT/complete_clusters.tsv")
    rows = csv.reader(clusters_csv, delimiter = '\t')
    next(rows)
    cmd_clusters = '''
    INSERT OR IGNORE INTO clusters (p100, p90, p30) VALUES (?,?,?)
    '''
    
    cmd_index_p90 = "CREATE INDEX clustersp90 ON clusters (p90)"
    cmd_index_p30 = "CREATE INDEX clustersp30 ON clusters (p30)"
    #cursor.execute(cmd7)
    cursor.execute("DROP TABLE CLUSTERS")
    cursor.execute("CREATE TABLE clusters (p100 text, p90 text, p30 text, PRIMARY KEY(p100));")
    cursor.executemany(cmd_clusters, rows)
    cursor.execute(cmd_index_p90)
    cursor.execute(cmd_index_p30)
    # cursor.execute(cmd1)
    # conn.execute(cmd2)
    # # cursor.execute(cmd3)
    # # cursor.execute(cmd4)
    # cursor.execute(cmd5)
    # cursor.execute(cmd6)
    conn.commit()
    conn.close()
move_clusters_to_genegraphdb()

