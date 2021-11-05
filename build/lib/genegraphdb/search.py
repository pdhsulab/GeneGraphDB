from collections import defaultdict
from genegraphdb import *
from neo4j import GraphDatabase
from Bio import SeqIO
import hashlib
import time
from os.path import abspath
import sqlite3
from cassandra.cluster import Cluster

#Global variables
hash_size = 25
fileName = "uniref100.fasta"
csv_results = "MinhashNeo4j_results_10_8.csv"
outfile_results = open(csv_results, "w")
#print('kmer_size,num_sequences,num_minhash,csv_time,minhash_time,protein_time,edges_time,total_time,search_time', file = outfile_results)

#Calculates hashcode given salt
def hashWithSalt(salt, kmer):
     return hashlib.sha256(str(salt).encode() + kmer.strip('*').encode()).hexdigest()[0:hash_size]
#Generates a set of kmers for a protein at the designated kmer size
def kmerSet(sequence, kmer_size):
    kmer_set = set()
    start_index = 0
    while(start_index + kmer_size <= len(sequence) + 1):
        kmerSeq = sequence[start_index : start_index + kmer_size]
        kmer_set.add(kmerSeq)
        start_index += 1
    return kmer_set
#Minhash of particular hash function for a set of kmers
def findMinHash(salt_num, kmers):
    min_hash = ""
    for kmer in kmers:
        test_min_hash = hashWithSalt(salt_num, kmer)
        if min_hash == "" or test_min_hash < min_hash:
            min_hash = test_min_hash
    return min_hash
#Writes protein to minhash pairs (num_min_hash possible pairs per protein)
def addSeqCSV(sequence, name, kmer_size, num_min_hash, outfile_phashes, outfile_minhashes, outfile_edges):
    kmers = kmerSet(sequence, kmer_size)
    print(name, len(sequence), sep = ",", file = outfile_phashes)
    for salt_num in range(0, num_min_hash):
        minhash = findMinHash(salt_num, kmers)
        #Writes to CSV
        print(minhash, sep = ",", file = outfile_minhashes)
        print(name, minhash, len(sequence), sep = ",", file = outfile_edges)


#Create reference Database for Neo4j w/ Minhashing
def createDB_Neo4j_MH(kmer_size, MAX_SEQUENCES, num_min_hash):
    #Opens and reads UniRef100 database
    tic_csv = time.time()
    fasta_sequences_uniref = SeqIO.parse(open(fileName), 'fasta')
    csv_path_edges = "minhash_gene_tmp.csv"
    outfile_edges = open(csv_path_edges, "w")
    csv_path_minhashes = "minhash_tmp.csv"
    outfile_minhashes = open(csv_path_minhashes, "w")
    csv_path_phashes = "protein_tmp.csv"
    outfile_phashes = open(csv_path_phashes, "w")
    print('name,minhash,plen', file = outfile_edges)
    print('name,plen', file = outfile_phashes)
    print('minhash', file = outfile_minhashes)
    num_sequences = 0
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        addSeqCSV(seq, hashcode, kmer_size, num_min_hash, outfile_phashes, outfile_minhashes, outfile_edges)
        if num_sequences >= MAX_SEQUENCES:
            break
        num_sequences += 1
    outfile_edges.close()
    outfile_phashes.close()
    outfile_minhashes.close()
    toc_csv = time.time()
    csv_time = toc_csv - tic_csv
    ##Neo4j Code
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    conn.query("CREATE CONSTRAINT minhash_id ON (n:Minhash) ASSERT n.hashid IS UNIQUE", db = "neo4j")
    #Clear previous data (remove if you're not testing)
    conn.query("MATCH (a) -[r] -> () DELETE a, r", db = "neo4j")
    conn.query("MATCH (a) delete a", db = "neo4j")
    #Load minhash nodes with timers
    tic_minhash = time.time()
    cmd_loadKmers = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      WITH DISTINCT row.minhash as minhash
      MERGE (m:Minhash {{hashid: minhash}})
      """.format(csv=abspath(csv_path_minhashes))
    conn.query(cmd_loadKmers, db = "neo4j")
    toc_minhash = time.time()
    minhash_time = toc_minhash - tic_minhash
    #Load protein nodes with timers
    tic_proteins = time.time()
    cmd_loadProteins = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      WITH DISTINCT row.name as pname, row.plen as plength
      MERGE (p:Protein {{name: pname, length: plength}})
      """.format(csv=abspath(csv_path_phashes))
    conn.query(cmd_loadProteins, db = "neo4j")
    toc_proteins = time.time()
    proteins_time = toc_proteins - tic_proteins
    #Load relationships with timers
    tic_edges = time.time()
    cmd_loadRelations = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      MATCH (m:Minhash {{hashid: row.minhash}}), (p:Protein {{name: row.name}})
      MERGE (m)-[r:MINHASH_OF]->(p)
      """.format(csv=abspath(csv_path_edges))
    conn.query(cmd_loadRelations, db = "neo4j")
    toc_edges = time.time()
    edges_time = toc_edges - tic_edges
    total_time = toc_edges - tic_csv

#Search reference database (Neo4j w/ Minhashing)
def searchDB_Neo4j_MH(kmer_input, num_min_hash, kmer_size):
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    #Sample kmer_input (for testing)
    kmer_input = "IMGPPDPILGMVRGVCIVEMGPPDPPSVTQFVVSFKAEFEDVPWLKEDIARSDAKQGSQLVALHLRLMAGGEVVNGLAKDLARQHGKSGALFYLNGLLNQSAPIFTQEPRAPYNTGNEDSPALLKIGTTIIENET"
    proteinToNum = defaultdict(lambda: 0)
    kmer_set = kmerSet(kmer_input, kmer_size)
    ticS = time.time()
    #Add all kmers for this protein to the Neo4j database
    for salt_num in range(num_min_hash):
        minhash = findMinHash(salt_num, kmer_set)
        cmd = """
            MATCH (a:Minhash) - [:MINHASH_OF] ->(p:Protein)
            WHERE a.hashid ='{query_minhash}'
            RETURN p.name
            """.format(query_minhash = minhash)
        matchingProteins = conn.query(cmd, db = "neo4j")
        if len( matchingProteins) != 0:
            for p in matchingProteins:
                proteinToNum[p] += 1
    tocS = time.time()
    search_time = tocS - ticS
    #Final results
    #print(kmer_size, MAX_SEQUENCES, num_min_hash, csv_time, minhash_time, proteins_time, edges_time, total_time, search_time, sep = ",", file = outfile_results)

#Add minhashes of sequence into SQL database
def addSeqMinHashes(hashcode, sequence, kmer_size, num_min_hash, c):
    kmer_set = kmerSet(sequence, kmer_size)
    tmp_insertions = 0
    for salt_num in range(0, num_min_hash):
        tmp_insertions += 1
        min_hash = findMinHash(salt_num, kmer_set)
        c.execute("""INSERT INTO minhashes
                     VALUES (?,?)""", (hashcode, min_hash))
    return tmp_insertions

#Create reference database for SQL w/ Minhashing
def createDB_SQL_MH(kmer_size, MAX_SEQUENCES, num_min_hash):
    #Variables
    rundate = "8_17_21"
    #hash_size = 25
    results_csv = "create_runtimes_" + rundate + ".csv"
    outfile_results = open(results_csv, "w")
    print('kmer_length,num_sequences,num_minhashes,insertion_time,index_time,total_time,num_rows,search_speed', file = outfile_results)

    #Create SQlite variables/table(s)
    conn = sqlite3.connect('SC_MH.db')
    c = conn.cursor()
    c.execute("""DROP TABLE minhashes""")
    c.execute("""CREATE TABLE IF NOT EXISTS minhashes
           (phash text,
            minhash text)""")
    #Variables
    referenceDB = "uniref100.fasta"
    num_sequences = 0
    num_insertions = 0
    #Add to SQL table
    tic = time.time()
    c.execute("PRAGMA synchronous = OFF")
    c.execute("PRAGMA journal_mode = MEMORY")
    c.execute("BEGIN TRANSACTION")
    fasta_sequences_uniref = SeqIO.parse(open(referenceDB), 'fasta')
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        num_insertions += addSeqMinHashes(hashcode, seq, kmer_size, num_min_hash, c)
        if num_sequences == MAX_SEQUENCES:
            break
        num_sequences += 1
    toc = time.time()
    c.execute("END TRANSACTION")
    insertion_time = toc - tic
    c.execute("BEGIN TRANSACTION")
    ticInd = time.time()
    c.execute("CREATE INDEX minhashID on minhashes(minhash)")
    tocInd = time.time()
    index_time = tocInd - ticInd
    total_time = tocInd - tic
    c.execute("END TRANSACTION")

def searchDB_SQL_MH(kmer_input, kmer_size, num_min_hash):
    conn = sqlite3.connect('SC_MH.db')
    c = conn.cursor()
    #Sample kmer_input (for testing)
    kmer_input = "IMGPPDPILGMVRGVCIVEMGPPDPPSVTQFVVSFKAEFEDVPWLKEDIARSDAKQGSQLVALHLRMAGGEVVNGLAKDLARQHGKSGALFYLNGLLNQSAPIFTQEPRAPYNTGNEDSPALLKIGTTIIENET"
    phash_to_num = defaultdict(lambda: 0)
    ticS = time.time()
    kmers = kmerSet(kmer_input, kmer_size)
    for salt_num in range(num_min_hash):
        minhash = findMinHash(salt_num, kmers)
        c.execute("SELECT phash FROM minhashes WHERE minhash = (?)", (minhash,))
        phashes = c.fetchall()
        for phash in phashes:
            phash_to_num[phash] += 1
    tocS = time.time()
    search_time = tocS - ticS
    #print(kmer_size, MAX_SEQUENCES, num_min_hash, insertion_time, index_time, total_time, num_insertions,search_time, sep =",", file = outfile_results)
    print("DONE")
    conn.close()

#Generates a set of kmers for a sequence at the designated kmer size
def kmerSet(sequence, kmer_size):
    kmer_set = set()
    start_index = 0
    while(start_index + kmer_size <= len(sequence) + 1):
        kmerSeq = sequence[start_index : start_index + kmer_size]
        kmer_set.add(kmerSeq)
        start_index += 1
    return kmer_set
#Add minhashes of sequence into SQL database
def addSeqKmers(hashcode, sequence, kmer_size, c):
    kmer_set = kmerSet(sequence, kmer_size)
    tmp_insertions = 0
    for kmer in kmer_set:
        tmp_insertions += 1
        c.execute("""INSERT INTO kmers
                     VALUES (?,?)""", (hashcode, kmer))
    return tmp_insertions

#SQL w/ Kmers
def createDB_SQL_Kmer(kmer_size, MAX_SEQUENCES):
    #"Global" variables
    referenceDB = "uniref100.fasta"
    csv_path = "kmerSQL_results_8_18.csv"
    outfile = open(csv_path, "w")
    #Add to SQL table
    conn = sqlite3.connect('SC.db')
    c = conn.cursor()
    c.execute("""DROP TABLE kmers""")
    c.execute("""CREATE TABLE IF NOT EXISTS kmers
                (phash text,
                 kmer text)""")
    c.execute("PRAGMA synchronous = OFF")
    c.execute("PRAGMA journal_mode = MEMORY")
    c.execute("BEGIN TRANSACTION")
    num_sequences = 0
    num_insertions = 0
    fasta_sequences_uniref = SeqIO.parse(open(referenceDB), 'fasta')
    tic = time.time()
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        num_insertions += addSeqKmers(hashcode, seq, kmer_size, c)
        if num_sequences == MAX_SEQUENCES:
            break
        num_sequences += 1
    c.execute("END TRANSACTION")
    toc = time.time()
    insertion_time = toc - tic
    #Create index on kmers (for searching after)
    c.execute("BEGIN TRANSACTION")
    ticInd = time.time()
    c.execute("CREATE INDEX kmerID on kmers(kmer)")
    c.execute("END TRANSACTION")
    tocInd= time.time()
    index_time = tocInd - ticInd
    total_time = tocInd - tic
def searchDB_SQL_Kmer(kmer_size):
    conn = sqlite3.connect('SC.db')
    c = conn.cursor()
    kmer_input = "IMGPPDPILGMVRGVCIVEMGPPDPPSVTQFVVSFKAEFEDVPWLKEDIARSDAKQGSQLVALHLRMAGGEVVNGLAKDLARQHGKSGALFYLNGLLNQSAPIFTQEPRAPYNTGNEDSPALLKIGTTIIENET"
    phash_to_num = defaultdict(lambda: 0)
    ticS = time.time()
    kmers = kmerSet(kmer_input, kmer_size)
    for kmer in kmers:
        c.execute("SELECT phash FROM kmers WHERE kmer = (?)", (kmer,))
        phashes = c.fetchall()
        for phash in phashes:
            phash_to_num[phash] += 1
    tocS = time.time()
    search_time = tocS - ticS
    #print(kmer_size, MAX_SEQUENCES, insertion_time, index_time, total_time, num_insertions, search_time, sep = ",", file = outfile)   


#Cassandra w/ Kmers
def createDB_Cass_Kmer(kmer_size, MAX_SEQUENCES):
    referenceDB = "uniref100.fasta"
    csv_path = "kmerCass_results_8_18.csv"
    outfile = open(csv_path, "w")
    cluster = Cluster()
    session = cluster.connect('kmers_test')
    session.execute("DROP TABLE kmersTable")
    createTableQry = """
    CREATE TABLE IF NOT EXISTS kmersTable(
         phash text,
         kmer text,
         primary key(kmer)
    );"""
    session.execute(createTableQry)
    num_sequences = 0
    num_insertions = 0
    fasta_sequences_uniref = SeqIO.parse(open(referenceDB), 'fasta')
    tic = time.time()
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        num_insertions += addSeqKmers(hashcode, seq, kmer_size, session)
        if num_sequences == MAX_SEQUENCES:
            break
        num_sequences += 1
    toc = time.time()
    insertion_time = toc - tic
def searchDB_Cass_Kmer(kmer_size):
    cluster = Cluster()
    session = cluster.connect('kmers_test')
    kmer_input = "IMGPPDPILGMVRGVCIVEMGPPDPPSVTQFVVSFKAEFEDVPWLKEDIARSDAKQGSQLVALHLRMAGGEVVNGLAKDLARQHGKSGALFYLNGLLNQSAPIFTQEPRAPYNTGNEDSPALLKIGTTIIENET"
    phash_to_num = defaultdict(lambda: 0)
    ticS = time.time()
    kmers = kmerSet(kmer_input, kmer_size)
    prepSearch = session.prepare("SELECT phash FROM kmersTable WHERE kmer = (?)")
    for kmer in kmers:
        phashes = session.execute(prepSearch, [kmer])
        for phash in phashes:
            phash_to_num[phash] += 1
    tocS = time.time()
    search_time = tocS - ticS
    #print(kmer_size, MAX_SEQUENCES, insertion_time, num_insertions, search_time, sep = ",", file = outfile)
