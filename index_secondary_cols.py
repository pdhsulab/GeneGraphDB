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


# In[1]:


def move_clusters_to_genegraphdb():
    conn = sqlite3.connect('genegraph.new.db')
    cursor = conn.cursor()
    cmd1 = "CREATE INDEX p2pwindow_p2hash ON prot2protwindow (p2hash)"
    cmd2 = "CREATE INDEX contig2sample_sampleid ON contig2sample (sampleid)"
    cmd3 = "CREATE INDEX crisprcoords_contighash ON crisprcoords (contighash)"
    cmd4 = "CREATE INDEX proteincoords_contighash ON proteincoords (contighash)"
    cmd5 = "CREATE INDEX p2p_p2hash ON prot2prot (p2hash)"
    cmd6 = "CREATE INDEX p2c_crisprhash ON prot2crispr (crisprhash)"
    cmd7 = "CREATE INDEX p2cwindow_crisprhash ON prot2crisprwindow (crisprhash)"
    cmd8 = "CREATE INDEX reppid ON stringent (reppid)"
    cmd9 = "CREATE INDEX reppid_permissive ON permissive (reppid)"
    cmd10 = "CREATE INDEX clustersp90 ON clusters (p90)"
    cmd11 = "CREATE INDEX clustersp30 ON clusters (p30)"
    cursor.execute(cmd10)
    cursor.execute(cmd11)
    cursor.execute(cmd1)
    cursor.execute(cmd2)
    cursor.execute(cmd3)
    cursor.execute(cmd4)
    cursor.execute(cmd5)
    cursor.execute(cmd6)
    cursor.execute(cmd7)
    cursor.execute(cmd8)
    cursor.execute(cmd9)

    conn.commit()
    conn.close()
move_clusters_to_genegraphdb()


# In[ ]:




