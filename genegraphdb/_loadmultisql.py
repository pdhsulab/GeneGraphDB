from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import proteinnode
from genegraphdb import proteinnodesql
from genegraphdb import crisprnode
from genegraphdb import _load
from genegraphdb import _loadsql
from genegraphdb import crisprnodesql
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import time

def _single(sample_id, google_bucket, distance, comment, outfile, samples_path=""):
    sample_id_path = sample_id + "/"
    fasta, protein, contigs = samples_path + sample_id_path + sample_id + ".fna.gz", \
                              samples_path + sample_id_path + sample_id + ".prodigal.faa.gz", \
                              samples_path + sample_id_path + sample_id + ".contigs.tsv.gz"
    tic = time.time()
    sorted_gff_name = crisprnodesql.merge_gff(sample_id, samples_path)
    _loadsql.load_proteins(sample_id, protein, samples_path)
    crisprid2crhash = crisprnodesql.load_CRISPRs(sample_id, samples_path)
    recid2contig = _loadsql.load_fasta(sample_id, fasta, contigs, samples_path)
    _loadsql.load_contig2sample(sample_id, contigs, samples_path)
    _loadsql.load_coords(sample_id, sorted_gff_name, recid2contig, crisprid2crhash, samples_path)
    proteinnodesql.create_all_protein_crispr_edge_csv(sample_id, distance,
                                                      samples_path)  # TO do - verify i need this
    toc = time.time()
    print("Loading the entire sql database, sans prot-prot and prot-crispr edges, took %f seconds" % (toc - tic) + "\n")
    if comment is None:
        comment = ""
    # print(sample_id + "," + str(toc - tic) + "," + comment)
    print(sample_id + "," + str(toc - tic) + ",null," + comment, file=outfile)

def bulk_connect_proteins_crisprs(max_distance, samples_path = ''):
    # to do - need to create protein2protein csvs from scratch
    print("Loading protein2protein edges...")
    tic = time.time()
    p2pcsv_paths, p2ccsv_paths, p2p_window_csv_paths, p2c_window_csv_paths= [], [], [], []
    for sample_id in os.listdir(samples_path):
        try:
            # to do - implement better way to check if the sample_id is actually a directory
            os.chdir(samples_path + sample_id)
            os.chdir("../..")
            sample_id_path = sample_id + "/"
            proteinnodesql.create_all_protein_crispr_edge_csv(sample_id, max_distance, samples_path)
            # to do: debug this here and in other loadmulti.py file
            p2pcsv_paths.append(samples_path + sample_id_path + "protein2protein.tmp.sql.csv")
            p2ccsv_paths.append(samples_path + sample_id_path + "protein2crispr.tmp.sql.csv")
            p2p_window_csv_paths.append(samples_path + sample_id_path + "protein2protein_window.tmp.sql.csv")
            p2c_window_csv_paths.append(samples_path + sample_id_path + "protein2crispr_window.tmp.sql.csv")
        except NotADirectoryError:
            print(sample_id + " is not a directory")
    os_cmd_p2p, os_cmd_p2c, os_cmd_p2p_window, os_cmd_p2c_window = ("awk 'FNR==1 && NR!=1{next;}{print}'",)*4
    # to do - start tic timer here. might be bottleneck to bulk loading
    for path in p2pcsv_paths:
        os_cmd_p2p += " " + path
    for path in p2ccsv_paths:
        os_cmd_p2c += " " + path
    for path in p2p_window_csv_paths:
        os_cmd_p2p_window += " " + path
    for path in p2c_window_csv_paths:
        os_cmd_p2c_window += " " + path

    os_cmd_p2p += " > " + samples_path + "bigfile_p2p.sql.csv"
    os_cmd_p2c += " > " + samples_path + "bigfile_p2c.sql.csv"
    os_cmd_p2p_window += " > " + samples_path + "bigfile_p2p_window.sql.csv"
    os_cmd_p2c_window += " > " + samples_path + "bigfile_p2c_window.sql.csv"

    os.system(os_cmd_p2p)
    os.system(os_cmd_p2c)
    os.system(os_cmd_p2p_window)
    os.system(os_cmd_p2c_window)
    # to do - stop timer here. print output
    proteinnodesql.load_csv(samples_path + "bigfile_p2p.sql.csv", samples_path + "bigfile_p2c.sql.csv",
                            samples_path + "bigfile_p2p_window.sql.csv", samples_path + "bigfile_p2c_window.sql.csv")
    toc = time.time()
    print("Loading protein2protein edges took %f seconds" % (toc - tic))
    return toc - tic