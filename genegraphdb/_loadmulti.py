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

def _single(sample_id, google_bucket, gene_neighbors, distance, comment, outfile):
    #outfile = outfile
    #os.chdir(sample_id)
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
    #proteinnode.connect_proteins_crisprs(distance, gene_neighbors)
    toc = time.time()
    #os.chdir("..")
    print("Loading the entire database took %f seconds" % (toc - tic) + "\n")
    if comment is None:
        comment = ""
    # print(sample_id + "," + str(toc - tic) + "," + comment)
    print(sample_id + "," + str(toc - tic) + ",null," + comment, file=outfile)

def bulk_load_protein_crispr_edges(distance, gene_neighbors):
    print("Loading protein2protein edges...")
    tic = time.time()
    protein2protein_csv = "protein2protein.tmp.csv"
    csv_paths = []
    for sample_id in os.listdir():
        try:
            # to do - implement better way to check if the sample_id is actually a directory
            os.chdir(sample_id)
            os.chdir("..")
            csv_paths.append(sample_id + "/protein2protein.tmp.csv")
            # os.system("cat file1.csv <(tail +2 file2.csv) <(tail +2 file3.csv) > bigfile.csv")
        except NotADirectoryError:
            print(sample_id + " is not a directory")
    os_cmd = "awk 'FNR==1 && NR!=1{next;}{print}'"
    for path in csv_paths:
        os_cmd += " " + path
    os_cmd += " > bigfile.csv"
    os.system(os_cmd)
    proteinnode.load_csv("bigfile.csv")
    toc = time.time()
    print("Loading protein2protein edges took %f seconds" % (toc - tic))
