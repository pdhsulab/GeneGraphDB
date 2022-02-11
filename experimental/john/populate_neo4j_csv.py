from db_inference import simple_neo4j
from utils import ggdb_logging
from utils import profile_util

# CLUSTERS_CSV = "file:///csv_exports/clusters_sample.csv"
# PROT2PROT_WINDOW_CSV = "file:///csv_exports/prot2protwindow_sample.csv"
CLUSTERS_CSV = "file:///csv_exports/clusters.csv"
PROT2PROT_WINDOW_CSV = "file:///csv_exports/prot2protwindow.csv"


def uniqueness_constraints(conn):
    conn.query("CREATE CONSTRAINT ON (p:P30) ASSERT p.p30 IS UNIQUE")
    conn.query("CREATE CONSTRAINT ON (p:P90) ASSERT p.p90 IS UNIQUE")
    conn.query("CREATE CONSTRAINT ON (p:P100) ASSERT p.p100 IS UNIQUE")


def load_clusters(conn):
    query = (
        """
            USING PERIODIC COMMIT 10000
            LOAD CSV WITH HEADERS FROM "%s" AS row
            MERGE (c30:P30 {p30: row.p30})
            MERGE (c90:P90 {p90: row.p90})
            MERGE (c100:P100 {p100: row.p100})
            MERGE (c30)-[:P30_CLUSTERING]->(c90)
            MERGE (c90)-[:P90_CLUSTERING]->(c100)
            RETURN count(*) as total
            """
        % CLUSTERS_CSV
    )
    conn.query(query)


def load_prot2prot(conn):
    query = (
        """
            USING PERIODIC COMMIT 10000
            LOAD CSV WITH HEADERS FROM "%s" AS row
            MERGE (n:P100 {p100: row.p1hash})
            MERGE (m:P100 {p100: row.p2hash})
            MERGE (n)-[:WINDOWED_NEIGHBOR]->(m)
            RETURN count(*) as total
            """
        % PROT2PROT_WINDOW_CSV
    )
    conn.query(query)


def main():
    conn = simple_neo4j.Neo4jConnection()
    with profile_util.time_monitor("Uniqueness constraints"):
        uniqueness_constraints(conn)

    with profile_util.time_monitor("Load clusters"):
        ggdb_logging.info("Loading clusters")
        load_clusters(conn)

    with profile_util.time_monitor("Load prot2protwindow"):
        ggdb_logging.info("Loading prot2protwindows")
        load_prot2prot(conn)

    ggdb_logging.info("Finished")


if __name__ == "__main__":
    main()
