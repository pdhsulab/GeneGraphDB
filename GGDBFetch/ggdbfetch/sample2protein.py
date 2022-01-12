from ggdbfetch import SAMPLE2PROTEIN_DB
from os.path import join
import sqlite3
from collections import defaultdict
from multiprocessing import Pool

def get_sample_to_p100s(p100s, dbpath, sample2path, threads):
    dbpath = join(dbpath, SAMPLE2PROTEIN_DB)

    args = [(p100, dbpath, sample2path) for p100 in p100s]
    with Pool(threads) as pool:
        results = pool.starmap(get_sample_id, args)
    sample2p100s = defaultdict(set)
    for res in results:
        for samp in res:
            for p100 in res[samp]:
                sample2p100s[samp].add(p100)

    return dict(sample2p100s)

def get_sample_id(p100, dbpath, sample2path):
    con = sqlite3.connect(dbpath)
    cur = con.cursor()
    cmd = 'SELECT sample_id FROM sample2protein WHERE p100="{}"'.format(p100)
    sample2p100s = defaultdict(set)
    for res in cur.execute(cmd):
        samp = res[0]
        if samp in sample2path:
            sample2p100s[res[0]].add(p100)
    cur.close()
    con.close()
    return sample2p100s
