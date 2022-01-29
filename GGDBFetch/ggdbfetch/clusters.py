import sqlite3
from os.path import join

from ggdbfetch import CLUSTERS_DB


def get_clusters(p30_id, dbpath):

    con = sqlite3.connect(join(dbpath, CLUSTERS_DB))
    cur = con.cursor()
    cmd = 'SELECT p100,p90 FROM clusters WHERE p30="{}"'.format(p30_id)

    p100_to_p90 = dict()
    p100_to_p30 = dict()
    p100s = set()
    for res in cur.execute(cmd):
        p100, p90 = res
        p100_to_p90[p100] = p90
        p100_to_p30[p100] = p30_id
        p100s.add(p100)

    cur.close()
    con.close()

    return p100s, p100_to_p90, p100_to_p30

