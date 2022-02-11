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
            ggdb_logging.warning("Query failed:", e)
        finally:
            if session is not None:
                session.close()
        return response
