from ggdbfetch import CLUSTERS_DB
from os.path import join
import sqlite3

def p30_to_p100(p30_id, dbpath):

    con = sqlite3.connect(join(dbpath, CLUSTERS_DB))
    cur = con.cursor()
    cmd = 'SELECT p100 FROM clusters WHERE p30="{}"'.format(p30_id)

    p100s = set()
    for res in cur.execute(cmd):
        p100s.add(res[0])

    cur.close()
    con.close()

    return p100s

