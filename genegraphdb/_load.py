from genegraphdb import *
from Bio import SeqIO
import gzip
import re
import sys
import os
from genegraphdb import graphdb

def _single(sample_id, fasta, protein, gff, contigs, google_bucket):

    load_proteins(protein)


def load_proteins(protein):

    outfile = open("proteins.tmp.csv", "w")
    print("ProteinId:ID,length:int,:LABEL", file=outfile)

    done = set()
    handle = gzip.open(protein, 'rt')
    for rec in SeqIO.parse(handle, 'fasta'):
        name_truncate = rec.id[:20]
        if name_truncate in done:
                continue
        done.add(name_truncate)
        length = len(str(rec.seq.strip('*')))
        print(name_truncate+","+str(length)+",Protein", file=outfile)
    outfile.close()

    cmd = "{neo4j_admin} import --database {db} --nodes {proteins}".format(
        neo4j_admin=NEO4J_ADMIN, db=DBNAME, proteins="proteins.tmp.csv"
    )
    os.system(cmd)

def load_fasta(fasta, contigs):
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)

    with gzip.open(contigs, 'rt') as infile:
        pass

    handle = gzip.open(protein, 'rt')
    for rec in SeqIO.parse(handle, 'fasta'):
        name_truncate = rec.id[:20]
        graphdb.add_node("Protein", {"hashid": name_truncate}, conn)

    handle.close()

