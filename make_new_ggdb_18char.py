#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import csv
import sqlite3
import time
from multiprocessing import Pool, cpu_count
import sys


# In[ ]:


def protein_stats_db():
    conn = sqlite3.connect('80kprotein_stats.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE '80kprotein_stats.old.db' AS db2"
    cmd2 = "INSERT OR IGNORE INTO proteins SELECT SUBSTRING(pid, 1, 18), isgapfree, length, iscomplete, sequence FROM db2.proteins"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
protein_stats_db()


# In[ ]:


def proteins():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO proteins SELECT SUBSTRING(hashid, 1, 18), length FROM db2.proteins"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#proteins()


# In[ ]:


def insert_rest():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO samples SELECT sampleid FROM db2.samples"
    cmd3 = "INSERT INTO crisprs SELECT SUBSTRING(hashid, 1, 18) FROM db2.crisprs"
    cmd4 = "INSERT INTO contigs SELECT SUBSTRING(hashid, 1, 18), length FROM db2.contigs"
    cmd5 = "INSERT INTO contig2sample SELECT SUBSTRING(contighashid, 1, 18), sampleid FROM db2.contig2sample"
    cmd6 = "INSERT INTO crisprcoords SELECT SUBSTRING(crisprhash, 1, 18), SUBSTRING(contighash, 1, 18), start, end FROM db2.crisprcoords"
    cmd7 = "INSERT INTO proteincoords SELECT SUBSTRING(phash, 1, 18), SUBSTRING(contighash, 1, 18), start, end, orientation FROM db2.proteincoords"
    cmd8 = "INSERT INTO prot2prot SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(p2hash, 1, 18) FROM db2.prot2prot"
    cmd9 = "INSERT INTO prot2crispr SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(crisprhash, 1, 18) FROM db2.prot2crispr"
    cmd10 = "INSERT INTO prot2protwindow SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(p2hash, 1, 18) FROM db2.prot2protwindow"
    cmd11 = "INSERT INTO prot2crisprwindow SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(crisprhash, 1, 18) FROM db2.prot2crisprwindow"
    cmd12 = "INSERT INTO stringent SELECT SUBSTRING(reppid, 1, 18), SUBSTRING(pid, 1, 18) FROM db2.stringent"
    cmd13 = "INSERT INTO permissive SELECT SUBSTRING(reppid, 1, 18), SUBSTRING(pid, 1, 18) FROM db2.permissive"
    cmd14 = "INSERT INTO clusters SELECT SUBSTRING(p100, 1, 18), SUBSTRING(p90, 1, 18), SUBSTRING(p30, 1, 18) FROM db2.clusters"
    
    conn.execute(cmd1)
    cursor.execute(cmd7)
    cursor.execute(cmd14)
    cursor.execute(cmd2)
    cursor.execute(cmd3)
    cursor.execute(cmd4)
    cursor.execute(cmd5)
    cursor.execute(cmd6)

    cursor.execute(cmd8)
    cursor.execute(cmd9)
    cursor.execute(cmd10)
    cursor.execute(cmd11)
    cursor.execute(cmd12)
    cursor.execute(cmd13)

    conn.commit()
    conn.close()

# insert_rest()


# In[ ]:


def samples():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO samples SELECT sampleid FROM db2.samples"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#samples()


# In[ ]:


def crisprs():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO crisprs SELECT SUBSTRING(hashid, 1, 18) FROM db2.crisprs"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#crisprs()


# In[ ]:


def contigs():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO contigs SELECT SUBSTRING(hashid, 1, 18), length FROM db2.contigs"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#contigs()


# In[ ]:


def contig2sample():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO contig2sample SELECT SUBSTRING(contighashid, 1, 18), sampleid FROM db2.contig2sample"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#contig2sample()


# In[ ]:


def crisprcoords():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO crisprcoords SELECT SUBSTRING(crisprhash, 1, 18), SUBSTRING(contighash, 1, 18), start, end FROM db2.crisprcoords"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#crisprcoords()


# In[ ]:


def proteincoords():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO proteincoords SELECT SUBSTRING(phash, 1, 18), SUBSTRING(contighash, 1, 18), start, end, orientation FROM db2.proteincoords"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#proteincoords()


# In[ ]:


def prot2prot():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO prot2prot SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(p2hash, 1, 18) FROM db2.prot2prot"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#prot2prot()


# In[ ]:


def prot2crispr():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO prot2crispr SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(crisprhash, 1, 18) FROM db2.prot2crispr"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#prot2crispr()


# In[ ]:


def prot2protwindow():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO prot2protwindow SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(p2hash, 1, 18) FROM db2.prot2protwindow"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#prot2protwindow()


# In[ ]:


def prot2crisprwindow():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO prot2crisprwindow SELECT SUBSTRING(p1hash, 1, 18), SUBSTRING(crisprhash, 1, 18) FROM db2.prot2crisprwindow"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#prot2crisprwindow()


# In[ ]:


def stringent():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO stringent SELECT SUBSTRING(reppid, 1, 18), SUBSTRING(pid, 1, 18) FROM db2.stringent"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#stringent()


# In[ ]:


def permissive():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO permissive SELECT SUBSTRING(reppid, 1, 18), SUBSTRING(pid, 1, 18) FROM db2.permissive"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#permissive()


# In[ ]:


def clusters():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "ATTACH DATABASE 'genegraph.db' AS db2"
    cmd2 = "INSERT INTO clusters SELECT SUBSTRING(p100, 1, 18), SUBSTRING(p90, 1, 18), SUBSTRING(p30, 1, 18) FROM db2.clusters"
    conn.execute(cmd1)
    cursor.execute(cmd2)
    conn.commit()
    conn.close()
#clusters()


# In[ ]:




