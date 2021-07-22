from genegraphdb import *
from neo4j import GraphDatabase
import time
import re
import hashlib
from os.path import abspath
from collections import defaultdict
import uuid

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



#Generates a set of kmers for a protein at the designated kmer size
def kmerSet(sequence, kmer_size):
    kmer_set = set()
    start_index = 0
    while(start_index + kmer_size <= len(sequence) + 1):
        kmerSeq = sequence[start_index : start_index + kmer_size]
        kmer_set.add(kmerSeq)
        start_index += 1
    return kmer_set

#Set up minhasing parameters/functions
hash_size = 25
num_min_hash = 25
num_to_salt = {}
for i in range(0,num_min_hash):
    num_to_salt[i] = uuid.uuid4().hex

#Calculates hashcode given salt
def hashWithSalt(salt, kmer):
    #Do this hashlib.pbkdf2_hmac('sha256', b'ACYAGYACYYG', b'1', 1).hex()
    return hashlib.sha256(salt.encode() + kmer.strip('*').encode()).hexdigest()[0:hash_size]

#Minhash of particular hash function for a set of kmers
def findMinHash(salt_num, kmers):
    min_hash = ""
    for kmer in kmers:
        test_min_hash = hashWithSalt(num_to_salt[salt_num], kmer)
        if min_hash == "" or test_min_hash < min_hash:
             min_hash = test_min_hash
    return min_hash

kmer_size = 5
def kmerdb():
    #Variable
    csv_path = "minhash_gene_tmp.csv"
    outfile = open(csv_path, "w")
    fileName = "uniref100.fasta"
    num_sequences = 0
    MAX_SEQUENCES = 1000

    #Writes protein to minhash pairs (num_min_hash possible pairs per protein)
    def addSeqCSV(sequence, name):
        kmers = kmerSet(sequence, kmer_size)
        for salt_num in range(0, num_min_hash):
            minhash = findMinHash(salt_num, kmers)
            #Writes to CSV
            print(name, minhash, len(sequence), sep = ",", file = outfile)
    
    #Opens and reads UniRef100 database
    tic_csv = time.time()
    print("STARTING CSV WRITING")
    fasta_sequences_uniref = SeqIO.parse(open(fileName), 'fasta')
    print('name,minhash,plen', file = outfile)
    for fasta in fasta_sequences_uniref:
        seq = str(fasta.seq)
        hashcode = hashlib.sha256(seq.strip('*').encode()).hexdigest()
        addSeqCSV(seq, hashcode)
        if num_sequences >= MAX_SEQUENCES:
            break
        if num_sequences % 10000 == 0:
            print(num_sequences)
        num_sequences += 1
    outfile.close()
    print("FINISHED WITH CSV")
    toc_csv = time.time()
    print("Writing CSV file took %f seconds" % (toc_csv - tic_csv))

    #Neo4j code
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)

    #Clear previous data 
    conn.query("MATCH (a) -[r] -> () DELETE a, r")
    conn.query("MATCH (a) delete a")

    #Load minhash nodes with timers
    tic_minhash = time.time()
    print("LOADING MINHASHES")
    cmd_loadKmers = """
    USING PERIODIC COMMIT
    LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
    WITH DISTINCT row.minhash as minhash
    MERGE (m:Minhash {{hashID: minhash}})
    """.format(csv=abspath(csv_path))
    conn.query(cmd_loadKmers, db = "neo4j")
    print("FINISHED LOADING MINHASHES")
    toc_minhash = time.time()
    print("Loading minhashes took %f seconds" % (toc_minhash-tic_minhash))

    #Load protein nodes with timers
    tic_proteins = time.time()
    print('LOADING PROTEINS')
    cmd_loadProteins = """
    USING PERIODIC COMMIT
    LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
    WITH DISTINCT row.name as pname, row.plen as plength
    MERGE (p:Protein {{name: pname, length: plength}})
    """.format(csv=abspath(csv_path))
    conn.query(cmd_loadProteins, db = "neo4j")
    print("FINISHED LOADING PROTEINS")
    toc_proteins = time.time()
    print("Loading proteins took %f seconds" % (toc_proteins - tic_proteins))

   #Load relationships with timers
    tic_edges = time.time()
    print("LOADING EDGES")
    cmd_loadRelations = """
    USING PERIODIC COMMIT
    LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
    MATCH (m:Minhash {{hashID: row.minhash}}), (p:Protein {{name: row.name}})
    MERGE (m)-[r:MINHASH_OF]->(p)
    """.format(csv=abspath(csv_path))
    conn.query(cmd_loadRelations, db = "neo4j")
    print('FINISHED LOADING EDGES')
    toc_edges = time.time()
    print("Loading edges took %f secondds" % (toc_edges - tic_edges))
    
    #Total number of nodes in Neo4j database (for testing)
    print(conn.query("MATCH (n) RETURN count(n)"))
    
    #Close connection
    conn.close()

#Takes in string of protein amino acid sequence
def search(protein):                                                                 
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    proteinToNum = defaultdict(lambda: 0)
    kmer_set = kmerSet(protein, kmer_size)

    #Add all kmers for this protein to the Neo4j database
    for salt_num in range(0, num_min_hash):
        min_hash = findMinHash(salt_num, kmer_set)
        cmd = """
           MATCH (m:Kmer) - [:MINHASH_OF] -> (p:Protein)
           WHERE m.hashID ='{minhash}'
           RETURN p.hashId
           """.format(minhash = min_hash)
        matchingProteins = conn.query(cmd, db = "neo4j")
        if matchingProteins is not None:
            for p in matchingProteins:
                proteinToNum[p] += 1
    print(proteinToNum)
    conn.close()



              
