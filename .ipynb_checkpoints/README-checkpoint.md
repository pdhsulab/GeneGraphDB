# GeneGraphDB
A command line tool to create and query a gene graph databases.
 
# Schema

## Nodes

### Sample
Types:
- Metagenome
  - A collection of contigs from many different organisms found in the same sample.
- Isolate
  - A collection of contigs belonging to a single organism that was isolated from a community.
- Metagenome-assembled genome (MAG)
  - A collection of contigs coming from a metagenome that are believed to belong to a single organism.

Abstract grouping of contigs.

### Contig
A contiguous sequence of nucleotides.

Properties:
contig_seq

### Protein
Amino acid translation of a predicted gene on a contig.

Properties:
protein_seq 

### StringentCluster

### PermissiveCluster

## Edges

### Contig2Sample

### Protein2Contig

### Protein2Protein

### Protein2Stringent

### Stringent2Permissive

