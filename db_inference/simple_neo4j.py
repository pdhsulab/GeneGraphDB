from typing import List

from neo4j import GraphDatabase

from utils import ggdb_logging

# NOTE: in developer_env/setup.sh we create a docker network, and both the neo4j container and usual ggdb container are
# run on this network. The shared network allows the containers to reference each other in URLs by container name.
# We specify the neo4j container name and authorization username/pw in setup.sh.
NEO4J_DOCKER_INSTANCE_NAME = "testneo4j"

DEFAULT_NEO4J_URI = f"bolt://{NEO4J_DOCKER_INSTANCE_NAME}:7687"
DEFAULT_NEO4J_USER = "neo4j"
DEFAULT_NEO4J_PWD = "test"

# Adapted from https://towardsdatascience.com/create-a-graph-database-in-neo4j-using-python-4172d40f89c4
class Neo4jConnection:
    def __init__(self, uri=DEFAULT_NEO4J_URI, user=DEFAULT_NEO4J_USER, pwd=DEFAULT_NEO4J_PWD):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            raise ConnectionError("Failed to create the driver:", e)

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
            # TODO: logging
            print("Query failed:", e)
        finally:
            if session is not None:
                session.close()
        return response


class SimpleNeo4j:
    def __init__(self):
        # TODO: validate data is as expected (e.g. that paper reader has loaded database)
        self.conn = Neo4jConnection()

    def get_p30_for_protein(self, p100_hash):
        query = f'MATCH (n:P30)-->(o:P90)-->(p:P100) WHERE p.p100 = "{p100_hash}" RETURN n.p30'
        resp = self.conn.query(query)
        p30_hash = resp[0][0]
        return p30_hash

    def get_num_p100s(self, p30_hash) -> int:
        query = f'MATCH (n:P30)-->(:P90)-->(p:P100) WHERE n.p30 = "{p30_hash}" RETURN count(p)'
        resp = self.conn.query(query)
        num_p100s = resp[0]["count(p)"]
        return num_p100s

    def get_targets_for_bait(self, p30_hash) -> List[str]:
        query = (
            "MATCH (bait:P30)-->(bn:P90)-->(bp:P100)--(tp:P100)<--(tn:P90)<--(tgt:P30) "
            f'WHERE bait.p30 = "{p30_hash}" RETURN collect(tgt.p30)'
        )
        resp = self.conn.query(query)
        tgts = resp[0][0]
        return tgts

    def get_num_shared(self, p30_A, p30_B):
        query = (
            "MATCH (n:P30)-->(o:P90)-->(p:P100)-[e]-(tp:P100)<--(to:P90)<--(tn:P30) "
            f'WHERE n.p30 = "{p30_A}" AND tn.p30="{p30_B}" RETURN count(e), count(p), count(tp)'
        )
        resp = self.conn.query(query)
        count_edges = resp[0]["count(e)"]
        count_prot_A = resp[0]["count(p)"]
        count_prot_B = resp[0]["count(tp)"]
        assert (
            count_edges == count_prot_A == count_prot_B
        ), f"Nonmatching counts {count_edges, count_prot_A, count_prot_B} for ({p30_A}, {p30_B})"
        return count_edges
