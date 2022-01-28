import re

from neo4j import GraphDatabase

from genegraphdb import *


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
    