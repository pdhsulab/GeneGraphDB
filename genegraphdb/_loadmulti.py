from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import proteinnode
from genegraphdb import crisprnode
from genegraphdb import _load
from Bio import SeqIO
import gzip
import re
import sys
import os
from os import remove
from os.path import abspath
import time

def _single(sample_id, google_bucket, distance, comment, outfile):
    sample_id_path = sample_id + "/"
    fasta, protein, contigs = sample_id_path + sample_id + ".fna.gz", \
                              sample_id_path + sample_id + ".prodigal.faa.gz", \
                              sample_id_path + sample_id + ".contigs.tsv.gz"
    tic = time.time()
    sorted_gff_name = crisprnode.merge_gff(sample_id)
    _load.load_proteins(sample_id, protein)
    crisprid2crhash = crisprnode.load_CRISPRs(sample_id)
    recid2contig = _load.load_fasta(sample_id, fasta, contigs)
    _load.load_contig2sample(sample_id, contigs)
    _load.load_coords(sample_id, sorted_gff_name, recid2contig, crisprid2crhash)
    toc = time.time()
    print("Loading the entire database, sans prot-prot and prot-crispr edges, took %f seconds" % (toc - tic) + "\n")
    if comment is None:
        comment = ""
    # print(sample_id + "," + str(toc - tic) + "," + comment)
    print(sample_id + "," + str(toc - tic) + ",null," + comment, file=outfile)

def bulk_connect_proteins_crisprs(max_distance):
    # to do - need to create protein2protein csvs from scratch
    print("Loading protein2protein edges...")
    tic = time.time()
    p2pcsv_paths, p2ccsv_paths, p2p_window_csv_paths, p2c_window_csv_paths= [], [], [], []
    for sample_id in os.listdir():
        try:
            # to do - implement better way to check if the sample_id is actually a directory
            os.chdir(sample_id)
            os.chdir("..")
            sample_id_path = sample_id + "/"
            # proteinnode.make_merged_coords_csv(sample_id)
            # outfile_prot_pair = open(sample_id_path + "protein2protein.tmp.csv", "w")
            # outfile_base_window = open(sample_id_path + "elemen2elemen_window.tmp.csv", "w")
            # outfile_prot_crispr_pair = open(sample_id_path + "protein2crispr.tmp.csv", "w")
            # print("recid,phash,qhash", file=outfile_prot_pair)
            # print("recid,phash,qhash", file=outfile_base_window)
            # print("recid,phash,qhash", file=outfile_prot_crispr_pair)
            #
            # merge_sorted_coords_csv = sample_id_path + "merged_sorted_coords.tmp.csv"
            # # write two distinct functions to create three types of protein edges (3 csvs)
            # proteinnode.create_protein_pair_csv(merge_sorted_coords_csv, outfile_prot_pair, outfile_prot_crispr_pair)
            # proteinnode.create_protein_window_csv(merge_sorted_coords_csv, max_distance, outfile_base_window)
            # outfile_prot_pair.close(), outfile_base_window.close(), outfile_prot_crispr_pair.close()
            proteinnode.create_all_protein_crispr_edge_csv(sample_id, max_distance)

            p2pcsv_paths.append(sample_id_path + "protein2protein.tmp.csv")
            p2ccsv_paths.append(sample_id_path + "protein2crispr.tmp.csv")
            p2p_window_csv_paths.append(sample_id_path + "protein2protein_window.tmp.csv")
            p2c_window_csv_paths.append(sample_id_path + "protein2crispr_window.tmp.csv")
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
    os_cmd_p2p += " > bigfile_p2p.csv"
    os_cmd_p2c += " > bigfile_p2c.csv"
    os_cmd_p2p_window += " > bigfile_p2p_window.csv"
    os_cmd_p2c_window += " > bigfile_p2c_window.csv"
    os.system(os_cmd_p2p)
    os.system(os_cmd_p2c)
    os.system(os_cmd_p2p_window)
    os.system(os_cmd_p2c_window)
    # to do - stop timer here. print output
    proteinnode.load_csv("bigfile_p2p.csv", "bigfile_p2c.csv", "bigfile_p2p_window.csv", "bigfile_p2c_window.csv")
    toc = time.time()
    print("Loading protein2protein edges took %f seconds" % (toc - tic))
    return toc - tic
