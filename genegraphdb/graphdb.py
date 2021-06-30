from genegraphdb import *
from neo4j import GraphDatabase
import re
import hashlib

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
    #ENTERPRISE ONLY: conn.query("CREATE CONSTRAINT exist_protein_hashid ON (n:Protein) ASSERT exists(n.hashid)", db=DBNAME)
    conn.query("CREATE CONSTRAINT uniq_contig_hashid ON (n:Contig) ASSERT n.hashid IS UNIQUE", db=DBNAME)
    #ENTERPRISE ONLY: conn.query("CREATE CONSTRAINT exist_contig_hashid ON (n:Contig) ASSERT exists(n.hashid)", db=DBNAME)
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

#Takes in string of protein amino acid sequence
def search(protein):
    start_index = 0                                                                  
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    kmer_size = 7
    proteinToNum = {}
    while(start_index + kmer_size <= len(protein)):
        kmerSeq = protein[start_index:start_index + kmer_size]
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
        start_index += 1
    print(proteinToNum)
    conn.close()

def kmerdb():
    #Variable
    csv_path = "kmer_protein_tmp.csv"
    outfile = open(csv_path, "w")
    fileName = "uniref100"
    kmer_size = 7
    num_sequences = 0
    MAX_SEQUENCES = 1000
    hash_to_sequence = {}
    
    #Add one sequence's kmers to CSV
    def addSeqCSV(sequence, hashcode):
        global num_insertions
        start_index = 0
        while (start_index + kmer_size < len(sequence)):
            kmer = sequence[start_index:start_index + kmer_size]
            ### Write to CSV
            print(kmer, hashcode, sep = ",", file = outfile)
            start_index += 1
            num_insertions += 1
    with open(fileName) as myFile:
        print('kmer,phash', file = outfile)
        for seq in myFile:
            #Calculate hashcode
            str = seq
            hashcode = hashlib.sha256(str.encode()).hexdigest()
            addSeqCSV(seq, hashcode)
            hash_to_sequence[hashcode] = seq
            if num_sequences >= MAX_SEQUENCES:
                break
            num_sequences += 1
    outfile.close()

    #Neo4j code
    conn = Neo4jConnection(DBURI, DBUSER, DBPASSWORD)
    #Clear previous data (used while testing)
    conn.query("MATCH (a) -[r] -> () DELETE a, r")
    conn.query("MATCH (a) delete a")
    #Load kmer nodes
    cmd_loadKmers = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      MERGE (k:Kmer {{kmerId: row.kmer}})
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadKmers, db = "neo4j")
    #Load protein nodes
    cmd_loadProteins = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      MERGE (p:Protein {{hashId: row.phash}})
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadProteins, db = "neo4j")
    #Load Relationships
    cmd_loadRelations = """
      USING PERIODIC COMMIT
      LOAD CSV WITH HEADERS FROM 'file:///{csv}' AS row
      MATCH (k:Kmer), (p:Protein)
      WHERE p.hashId = row.phash AND k.kmerId = row.kmer
      MERGE (k)-[r:KMER_OF]->(p
      """.format(csv=abspath(csv_path))
    conn.query(cmd_loadRelations, db = "neo4j")

    conn.close()


              
