from genegraphdb import *
from neo4j import GraphDatabase
import time
import re
import hashlib
from os.path import abspath
from collections import defaultdict

class Neo4jConnection:

    def __init__(self, uri, user, pwd):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            print("Failed to create the driver:", e)

    def close(self):
        if self.__driver is not None:
            self.__driver.close()

    def query(self, query, parameters=None, db=None):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        try:
            session = self.__driver.session(database=db) if db is not None else self.__driver.session()
            response = list(session.run(query, parameters))
        except Exception as e:
            print("Query failed:", e)
        finally:
            if session is not None:
                session.close()
        return response

def showdatabases():
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    return conn.query("SHOW DATABASES")
    conn.close()

def hasdb():
    for rec in showdatabases():
        if rec['name'] == DBNAME:
            return True
    return False

def createdb():
    print("Creating database: %s" % DBNAME)

    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    conn.query("CREATE DATABASE %s" % DBNAME)
    conn.query("CREATE CONSTRAINT uniq_protein_hashid ON (n:Protein) ASSERT n.hashid IS UNIQUE", db=DBNAME)
    # ENTERPRISE ONLY: conn.query("CREATE CONSTRAINT exist_protein_hashid ON (n:Protein) ASSERT exists(n.hashid)", db=DBNAME)
    conn.query("CREATE CONSTRAINT uniq_contig_hashid ON (n:Contig) ASSERT n.hashid IS UNIQUE", db=DBNAME)
     # ENTERPRISE ONLY: conn.query("CREATE CONSTRAINT exist_contig_hashid ON (n:Contig) ASSERT exists(n.hashid)", db=DBNAME)
    conn.close()

    print("Database created")

def num_nodes():
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    num_nodes = conn.query("MATCH (n) RETURN count(*)", db=DBNAME)[0]['count(*)']
    conn.close()

    return num_nodes

def num_rels():
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    num_rels = conn.query("MATCH ()-[r]->() RETURN count(r) as count", db=DBNAME)[0]['count']
    conn.close()

    return num_rels


def add_node(nodetype, properties=None, conn=None):

    close = False
    if conn is None:
        conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
        close = True

    if properties is not None:
        cmd = "CREATE (n:{ntype} {props})".format(ntype=nodetype, props=str(properties))
        cmd = re.sub(r"'([^':]+)':", r'\1:', cmd)
    else:
        cmd = "CREATE (n:{ntype} )".format(ntype=nodetype)

    conn.query(cmd, db=DBNAME)

    if close:
        conn.close()

def add_nodes(nodes, conn=None):
    close = False
    if conn is None:
        conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
        close = True

    cmd = "MERGE " + ",".join(nodes)
    conn.query(cmd, db=DBNAME)

    if close:
        conn.close()
    pass

def protein_exists(hashid, conn=None):

    close = False
    if conn is None:
        conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
        close = True

    cmd = "MATCH (n:Protein {hashid: '%s'}) return count(*) as count" % hashid
    exists = conn.query(cmd, db=DBNAME)[0]['count'] > 0

    if close:
        conn.close()

    return exists

#Generates a set of kmers for the desired protein at the desired kmer size
def kmerSet(protein, kmer_size):
    kmer_set = set()
    start_index = 0
    while(start_index + kmer_size <= len(protein) + 1):
        kmerSeq = protein[start_index:start_index + kmer_size]
        kmer_set.add(kmerSeq)
        start_index += 1
    return kmer_set

#Takes in string of protein amino acid sequence
def search(protein):                                                                 
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    kmer_size = 7
    proteinToNum = {}
    kmer_set = kmerSet(protein, kmer_size)

    #Add all kmers for this protein to the Neo4j database
    for kmerSeq in kmer_set:
        print(kmerSeq)
        cmd = """
           MATCH (a:Kmer) - [:KMER_OF] ->(p:Protein)
           WHERE a.kmerId ='{kmerName}'
           RETURN p.hashId
           """.format(kmerName=kmerSeq)
        matchingProteins = conn.query(cmd, db = "neo4j")
        if matchingProteins is not None:
            for p in matchingProteins:
                if p in proteinToNum:
                    proteinToNum[p] += 1
                else:
                    proteinToNum[p] = 1
    print(proteinToNum)
    conn.close()

#Add one sequence's kmers to CSV
def addSeqCSV(sequence, hashcode, kmer_size, outfile):
    kmer_set = kmerSet(sequence, kmer_size)
    for kmer in kmer_set:
        ### Write to CSV
        print(kmer, hashcode, sep = ",", file = outfile)

def kmerdb():
    #Variable
    csv_path = "kmer_protein_tmp.csv"
    outfile = open(csv_path, "w")
    fileName = "uniref100"
    kmer_size = 7
    num_sequences = 0
    MAX_SEQUENCES = 1000
    
    tic_csv = time.time()
    print('STARTING CSV WRITING')
    fasta_sequences_uniref = SeqIO.parse(open(fileName), 'fasta')
    print('kmer,phash', file = outfile)
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        addSeqCSV(seq, hashcode, kmer_size, outfile)
        if num_sequences >= MAX_SEQUENCES:
            break
        if num_sequences%100 == 0:
            print('FINISHED WITH {0} SEQUENCES'.format(num_sequences))
        num_sequences += 1
    outfile.close()
    print('FINISHED WITH CSV')
    toc_csv = time.time()
    print("Writing CSV file took %f seconds" % (toc_csv - tic_csv))

    #Neo4j code
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    #Clear previous data (used while testing)
    conn.query("MATCH (a) -[r] -> () DELETE a, r")
    conn.query("MATCH (a) delete a")
    #Load kmer nodes
    tic_kmers = time.time()
    print("LOADING KMERS")
    cmd_loadKmers = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      WITH DISTINCT row.kmer as kmer
      MERGE (k:Kmer {{kmerId: kmer}})
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadKmers, db = "neo4j")
    print('FINISHED LOADING KMERS')
    toc_kmers = time.time()
    print("Loading kmers took %f seconds" % (toc_kmers-tic_kmers))
    #Load protein nodes
    tic_proteins = time.time()
    print('LOADING PROTEINS')
    cmd_loadProteins = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      WITH DISTINCT row.phash as phash
      MERGE (p:Protein {{hashId: phash}})
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadProteins, db = "neo4j")
    print("FINISHED LOADING PROTEINS")
    toc_proteins = time.time()
    print("Loading proteins took %f seconds" % (toc_proteins - tic_proteins))
    #Load Relationships
    tic_edges = time.time()
    print('LOADING EDGES')
    cmd_loadRelations = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      MATCH (k:Kmer), (p:Protein)
      WHERE p.hashId = row.phash AND k.kmerId = row.kmer
      MERGE (k)-[r:KMER_OF]->(p
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadRelations, db = "neo4j")
    print('FINISHED LOADING EDGES')
    toc_edges = time.time()
    print("Loading edges took %f secondds" % (toc_edges - tic_edges))

    print(conn.query("MATCH (n) RETURN count(n)"))
    conn.close()


              
