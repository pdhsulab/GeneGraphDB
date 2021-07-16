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
    outfile = outfile
    os.chdir(sample_id)
    fasta, protein, contigs = sample_id + ".fna.gz", sample_id + ".prodigal.faa.gz", sample_id + ".contigs.tsv.gz"
    tic = time.time()
    sorted_gff_name = crisprnode.merge_gff(sample_id)
    _load.load_proteins(protein)
    crisprid2crhash = crisprnode.load_CRISPRs(sample_id)
    recid2contig = _load.load_fasta(fasta, contigs)
    _load.load_contig2sample(sample_id, contigs)
    _load.load_coords(sorted_gff_name, recid2contig, crisprid2crhash)
    #proteinnode.connect_proteins_crisprs(distance, gene_neighbors)
    toc = time.time()
    os.chdir("..")
    print("single sample id " + sample_id)
    print("Loading the entire database took %f seconds" % (toc - tic) + "\n")
    print(sample_id + "," + str(toc - tic) + "," + comment)
    print(sample_id + "," + str(toc - tic) + "," + comment, file=outfile)

def bulk_load_protein_crispr_edges():
    print("bulk edges loaded")


