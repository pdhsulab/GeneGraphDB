from ggdbfetch import CLUSTERS_DB
from os.path import join
import sqlite3
from collections import defaultdict

def get_sample_to_p100s(p100s, dbpath):
    con = sqlite3.connect(join(dbpath, CLUSTERS_DB))
    cur = con.cursor()

    sample2p100s = defaultdict(set)
    for p100 in p100s:
        cmd = 'SELECT sample_id FROM sample2protein WHERE p100="{}"'.format(p100)
        for res in cur.execute(cmd):
            sample2p100s[res[0]].add(p100)

    cur.close()
    con.close()

    return sample2p100s