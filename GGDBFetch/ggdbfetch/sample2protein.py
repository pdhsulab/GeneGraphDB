import sqlite3
from collections import defaultdict
from os.path import join

from ggdbfetch import SAMPLE2PROTEIN_DB
from tqdm import tqdm


def get_sample_to_p100s(p100s, dbpath):
    con = sqlite3.connect(join(dbpath, SAMPLE2PROTEIN_DB))
    cur = con.cursor()

    sample2p100s = defaultdict(set)
    for p100 in tqdm(p100s):
        cmd = 'SELECT sample_id FROM sample2protein WHERE p100="{}"'.format(p100)
        for res in cur.execute(cmd):
            samp = res[0]
            sample2p100s[res[0]].add(p100)
    cur.close()
    con.close()

    return dict(sample2p100s)

