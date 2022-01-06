# GeneGraphDB
A command line tool to create and query a gene graph databases.
 
# Schema

## Nodes

### samples
Types:
- Metagenome
  - A collection of contigs from many different organisms found in the same sample.
- Isolate
  - A collection of contigs belonging to a single organism that was isolated from a community.
- Metagenome-assembled genome (MAG)
  - A collection of contigs coming from a metagenome that are believed to belong to a single organism.

Abstract grouping of contigs.

CREATE TABLE samples (sampleid text, PRIMARY KEY (sampleid));

### contigs
A contiguous sequence of nucleotides.

CREATE TABLE contigs (hashid text, length real, PRIMARY KEY (hashid, length));

### proteins
Amino acid translation of a predicted gene on a contig.

CREATE TABLE proteins (hashid text, length real, PRIMARY KEY(hashid, length));

### crisprs
DNA sequence of CRISPR arrays

CREATE TABLE crisprs (hashid text, PRIMARY KEY (hashid));

### clusters
CREATE TABLE stringent (reppid text, pid text, PRIMARY KEY(pid));
- pid is p100 protein, reppid is p90 cluster representative

CREATE TABLE permissive (reppid text, pid text, PRIMARY KEY(pid));
- pid is p90 cluster representative, reppid is p30 cluster representative

## Edges

### Contig2Sample
CREATE TABLE contig2sample (contighashid text, sampleid text, PRIMARY KEY (contighashid, sampleid));

### Protein2Contig, Crispr2Contig
CREATE TABLE crisprcoords (crisprhash text, contighash text, start real, end real, PRIMARY KEY (crisprhash, contighash, start));

CREATE TABLE proteincoords (phash text, contighash text, start real, end real, orientation text, PRIMARY KEY (phash, contighash, start));

### Protein2Protein, Protein2Crispr
CREATE TABLE prot2prot (p1hash text, p2hash text, PRIMARY KEY (p1hash, p2hash));

CREATE TABLE prot2crispr (p1hash text, crisprhash text, PRIMARY KEY (p1hash, crisprhash));

CREATE TABLE prot2protwindow (p1hash text, p2hash text, PRIMARY KEY (p1hash, p2hash));
- window size is 5kb

CREATE TABLE prot2crisprwindow (p1hash text, crisprhash text, PRIMARY KEY (p1hash, crisprhash));
- window size is 5kb

### Protein2Stringent2permissive
CREATE TABLE clusters (p100, p90, p30, PRIMARY KEY(p100));

# Constructing GGDB 

## Making most nodes / relationships
Run genegraphdb package:

ggdb load multisql ...

## Making cluster nodes / relationships
### Run scripts in this order
get_mmseqs_input.ipynb
- Make multifasta files for mmseqs input

01.stringent_clusters.sh 
- Make stringent clusters

get_testdb_protein_stats.ipynb
- Make 80kprotein_stats.db

get_stringent_cluster_rep.ipynb
- Define which proteins are good stringent clusters
- Output clu_rep_stringent_final.csv
- Output clu_badrep_stringent.csv

make_permissive_clusters.ipynb
- Output clu_perm_mmseqs_input.faa
- This multifasta file is input for mmseqs

02.permissive_clusters.sh
- Make permissive clusters

03.combine_clusters.py
- Combines stringent and permissive clusters into one .csv file, complete_clusters.tsv
