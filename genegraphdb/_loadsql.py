from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import proteinnode
from genegraphdb import proteinnodesql
from genegraphdb import crisprnode
from genegraphdb import crisprnodesql
from genegraphdb import testing
from Bio import SeqIO
import gzip
import re
import sys
import csv
import os
from os import remove
from os.path import abspath
import time
import sqlite3

def _single(sample_id, google_bucket, distance, comment, outfilename, samples_path = '', clean_files=False):
    outfile = open(outfilename, "a") # to do - is this the only way to make multithreading work? best practices?
    try:
        sample_id_path = sample_id + "/"
        fasta, protein, contigs = samples_path + sample_id_path + sample_id + ".fna.gz", \
                                  samples_path + sample_id_path + sample_id + ".prodigal.faa.gz", \
                                  samples_path + sample_id_path + sample_id + ".contigs.tsv.gz"
        tic = time.time()
        sorted_gff_name = crisprnodesql.merge_gff(sample_id, samples_path)
        load_proteins(sample_id, protein, samples_path)
        crisprid2crhash = crisprnodesql.load_CRISPRs(sample_id, samples_path)
        recid2contig = load_fasta(sample_id, fasta, contigs, samples_path)
        load_contig2sample(sample_id, contigs, samples_path)
        load_coords(sample_id, sorted_gff_name, recid2contig, crisprid2crhash, samples_path)
        p2p_edge_load_time = proteinnodesql.connect_proteins_crisprs(sample_id, distance, samples_path)
        toc = time.time()
        print("Loading the entire database took %f seconds" % (toc - tic) + "\n")
        if comment is None:
            comment = ""
        if clean_files:
            testing.clean_files(sample_id, samples_path)
        print(sample_id + "," + str(toc - tic) + "," + str(p2p_edge_load_time) + "," + comment, file=outfile)
    except Exception as e:
        testing.log_errors_multisql_loadsql(samples_path, sample_id, e)
    outfile.close()

def load_proteins(sample_id, protein, samples_path = ''):
    print("Loading proteins...")
    tic = time.time()
    outfile = open(samples_path + sample_id + "/proteins.tmp.sql.csv", "w")
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

    con = sqlite3.connect('genegraph.db')
    cur = con.cursor()
    protein_csv_path = samples_path + sample_id + '/proteins.tmp.sql.csv'
    protein_csv = open(protein_csv_path)
    rows = csv.reader(protein_csv)
    next(rows)
    cmd = '''
    INSERT OR IGNORE INTO proteins (hashid, length) VALUES (?,?)
    '''
    cur.executemany(cmd, rows)
    con.commit()
    con.close()

    toc = time.time()
    #remove(samples_path + sample_id + '/proteins.tmp.sql.csv')
    print("Loading proteins took %f seconds" % (toc-tic))

def load_fasta(sample_id, fasta, contigs, samples_path = ''):
    print("Loading contigs...")
    tic = time.time()
    recid2contig = dict()
    with gzip.open(contigs, 'rt') as infile:
        infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            recid2contig[line[0]] = (line[1], int(line[3]))

    outfile = open(samples_path + sample_id + "/contigs.tmp.sql.csv", "w")
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

    con = sqlite3.connect('genegraph.db')
    cur = con.cursor()
    contig_csv_path = samples_path + sample_id + '/contigs.tmp.sql.csv'
    contig_csv = open(contig_csv_path)
    rows = csv.reader(contig_csv)
    next(rows)
    #cur.execute("PRAGMA journal_mode=WAL")
    cmd = '''
    INSERT OR IGNORE INTO contigs (hashid, length) VALUES (?,?)
    '''
    cur.executemany(cmd, rows)
    con.commit()
    con.close()
    toc = time.time()

    print("Loading contigs took %f seconds" % (toc-tic))
    # remove(samples_path + sample_id + '/contigs.tmp.sql.csv')

    return recid2contig

def load_contig2sample(sample_id, contigs, samples_path = ''):
    print("Loading contig2sample edges...")
    tic = time.time()
    outfile = open(samples_path + sample_id + "/contig2sample.tmp.sql.csv", "w")
    print("chash,sampleid", file=outfile)
    with gzip.open(contigs, 'rt') as infile:
        infile.readline()
        for line in infile:
            line = line.strip().split('\t')
            print(line[1][:20]+','+str(sample_id), file=outfile)
    outfile.close()
    con = sqlite3.connect('genegraph.db')
    cur = con.cursor()
    contig2sample_csv_path = samples_path + sample_id + '/contig2sample.tmp.sql.csv'
    contig2sample_csv = open(contig2sample_csv_path)
    rows = csv.reader(contig2sample_csv)
    next(rows)
    #cur.execute("PRAGMA journal_mode=WAL")
    cmd = '''
    INSERT OR IGNORE INTO contig2sample (contighashid, sampleid) VALUES (?,?)
    '''
    cur.executemany(cmd, rows)
    cmd2 = '''
    INSERT OR IGNORE INTO samples (sampleid) VALUES ("{}")
    '''.format(sample_id)
    cur.execute(cmd2)

    con.commit()
    con.close()
    toc = time.time()

    print("Loading contig2sample edges took %f seconds" % (toc-tic))

    #remove(samples_path + sample_id + "/contig2sample.tmp.csv")

def load_coords(sample_id, gff, recid2contig, crisprid2crhash, samples_path = ''):
    print("Loading gene and crispr coords...")
    tic = time.time()
    with open(samples_path + sample_id + "/gene_coords.tmp.sql.csv", "w") as g_out, \
            open(samples_path + sample_id + "/crispr_coords.tmp.sql.csv", "w") as cr_out:
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
    crisprnodesql.load_crispr_coords(sample_id, samples_path)
    # To do - debug this. why does gene_coords have a repeat? How can I get SQL to skip vs error out
    proteinnodesql.load_protein_coords(sample_id, samples_path)
    toc = time.time()

    print("Loading gene and crispr coords took %f seconds" % (toc-tic))