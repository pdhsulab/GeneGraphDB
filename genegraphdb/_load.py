from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import proteinnode
from genegraphdb import crisprnode
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import time

def _single(sample_id, google_bucket, gene_neighbors, distance, comment):
    outfile = open("ggdb_load_stats.csv", "a")
    print("sample_id,load_time,comment")
    os.chdir(sample_id)
    fasta, protein, contigs = sample_id + ".fna.gz", sample_id + ".prodigal.faa.gz", sample_id + ".contigs.tsv.gz"
    tic = time.time()
    sorted_gff_name = crisprnode.merge_gff(sample_id)
    load_proteins(protein)
    crisprid2crhash = crisprnode.load_CRISPRs(sample_id)
    recid2contig = load_fasta(fasta, contigs)
    load_contig2sample(sample_id, contigs)
    load_coords(sorted_gff_name, recid2contig, crisprid2crhash)
    proteinnode.connect_proteins_crisprs(distance, gene_neighbors)
    toc = time.time()
    os.chdir("..")
    print("Loading the entire database took %f seconds" % (toc - tic) + "\n")
    print(sample_id + "," + str(toc - tic) + "," + comment, file=outfile)


def load_proteins(protein):
    print("Loading proteins...")
    tic = time.time()
    outfile = open("proteins.tmp.csv", "w")
    print("hashid,length", file=outfile)

    done = set()
    handle = gzip.open(protein, 'rt')
    for rec in SeqIO.parse(handle, 'fasta'):
        name_truncate = rec.id[:20]
        if name_truncate in done:
                continue
        done.add(name_truncate)
        length = len(str(rec.seq.strip('*')))
        print(name_truncate+","+str(length), file=outfile)
    outfile.close()
    del done

    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = "LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row MERGE (n:Protein {{hashid: row.hashid, length:row.length}})".format(
        csv=abspath('proteins.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()
    toc = time.time()
    remove('proteins.tmp.csv')
    print("Loading proteins took %f seconds" % (toc-tic))



def load_fasta(fasta, contigs):
    print("Loading contigs...")
    tic = time.time()
    recid2contig = dict()
    with gzip.open(contigs, 'rt') as infile:
        infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            recid2contig[line[0]] = (line[1], int(line[3]))

    outfile = open("contigs.tmp.csv", "w")
    print("hashid,length", file=outfile)

    done = set()
    handle = gzip.open(fasta, 'rt')
    for rec in SeqIO.parse(handle, 'fasta'):
        contig, clength = recid2contig[rec.id]
        name_truncate = contig[:20]
        if name_truncate in done:
            continue
        done.add(name_truncate)
        length = len(str(rec.seq.strip('*')))
        print(name_truncate+","+str(clength), file=outfile)
    outfile.close()
    del done

    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd = """LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
          MERGE (n:Contig {{hashid: row.hashid, length:row.length}})
          """.format(
        csv=abspath('contigs.tmp.csv')
    )
    print(cmd)
    conn.query(cmd, db=DBNAME)
    conn.close()
    toc = time.time()

    print("Loading contigs took %f seconds" % (toc-tic))
    remove('contigs.tmp.csv')

    return recid2contig

def load_contig2sample(sample_id, contigs):
    print("Loading contig2sample edges...")
    tic = time.time()
    outfile = open("contig2sample.tmp.csv", "w")
    print("chash,sampleid", file=outfile)
    with gzip.open(contigs, 'rt') as infile:
        infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            print(line[1][:20]+','+str(sample_id), file=outfile)
    outfile.close()
    conn = graphdb.Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    cmd_make_node = """
                   MERGE (s:Sample {{sampleID: "{id}"}})
                    """.format(id=str(sample_id))
    cmd_load_edges = """
          LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row 
          MATCH (c:Contig), (s:Sample)
          WHERE c.hashid = row.chash AND s.sampleID = row.sampleid
          MERGE (c)-[e:contig2sample]->(s)
          """.format(
        csv=abspath('contig2sample.tmp.csv')
    )
    print(cmd_make_node)
    conn.query(cmd_make_node, db=DBNAME)
    print(cmd_load_edges)
    conn.query(cmd_load_edges, db=DBNAME)
    conn.close()
    toc = time.time()

    print("Loading contig2sample edges took %f seconds" % (toc-tic))

    #remove("contig2sample.tmp.csv")

def load_coords(gff, recid2contig, crisprid2crhash):
    print("Loading gene and crispr coords...")
    tic = time.time()
    with open("gene_coords.tmp.csv", "w") as g_out, open("crispr_coords.tmp.csv", "w") as cr_out:
        print("recid,phash,chash,start,end,orient", file=g_out)
        print("recid,crhash,chash,start,end", file=cr_out)
        with open(gff, "rt") as infile:
            for line in infile:
                if line.startswith("#"):
                    continue
                line = line.strip().split('\t')
                if len(line) != 9:
                    print(len(line))
                    continue
                if "Prodigal" in line[1]:
                    phash = line[-1].split("ID=")[-1].split(';')[0][:20]
                    chash = recid2contig[line[0]][0][:20]
                    print(line[0], phash, chash, line[3], line[4], line[6], sep=',', file=g_out)
                elif "minced" in line[1]:
                    old_id = line[-1].split("ID=")[-1].split(';')[0][:20]
                    crhash = crisprid2crhash[old_id]
                    chash = recid2contig[line[0]][0][:20]
                    print(line[0], crhash, chash, line[3], line[4], sep=',', file=cr_out)
    crisprnode.load_crispr_coords()
    proteinnode.load_protein_coords()
    toc = time.time()

    print("Loading gene and crispr coords took %f seconds" % (toc-tic))

    #remove("gene_coords.tmp.csv")
    #remove("crispr_coords.tmp.csv")