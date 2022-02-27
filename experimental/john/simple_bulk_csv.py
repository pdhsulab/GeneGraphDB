import os

from db_inference import simple_sql_db
from utils import ggdb_logging, profile_util

CLUSTER_LIMIT = 1000000  # max is ~200M
NEIGHBOR_LIMIT = 1000000  # max is ~1B


def write_neo4j_index_to_file(neo4j_set, fpath, node_name, property_name):
    with open(fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing {fpath}"):
        csv_fp.write(f"{property_name}:ID({node_name})\n")
        for protein_hash in neo4j_set:
            csv_fp.write(f"{protein_hash}\n")


def cluster_edges_to_file(px_to_py_set, px_node_name, py_node_name, fpath):
    with open(fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing {fpath}"):
        csv_fp.write(f":START_ID({px_node_name}),:END_ID({py_node_name})\n")
        for px_py_id in px_to_py_set:
            px_id, py_id = px_py_id.split("_")
            csv_fp.write(f"{px_id},{py_id}\n")


def main():
    ggdb_logging.info(f"CLUSTER_LIMIT is {CLUSTER_LIMIT}")
    ggdb_logging.info(f"NEIGHBOR_LIMIT is {NEIGHBOR_LIMIT}")

    cluster_query = "SELECT * FROM clusters"
    if CLUSTER_LIMIT is not None:
        cluster_query += f" LIMIT {CLUSTER_LIMIT}"

    neighbor_query = "SELECT * FROM prot2protwindow"
    if NEIGHBOR_LIMIT is not None:
        neighbor_query += f" LIMIT {NEIGHBOR_LIMIT}"

    sql_db = simple_sql_db.SimpleSqlDb()
    csv_dir = "/GeneGraphDB/data/neo4j/import"
    p100_fpath = os.path.join(csv_dir, "p100.csv")
    p100_neighbor_fpath = os.path.join(csv_dir, "p100_prot2prot_window_p100.csv")
    p90_fpath = os.path.join(csv_dir, "p90.csv")
    p90_p100_fpath = os.path.join(csv_dir, "p90_cluster_p100.csv")
    p30_fpath = os.path.join(csv_dir, "p30.csv")
    p30_p90_fpath = os.path.join(csv_dir, "p30_cluster_p90.csv")

    unique_p30s = set()
    unique_p90s = set()
    unique_p100s = set()

    # read through clusters file to gather unique p30s, p90s, p100s
    ggdb_logging.info("Reading clusters for unique proteins")
    num_read_cluster = 0
    cur = sql_db.conn.cursor()
    with profile_util.memory_monitor("clusters_unique_proteins"), profile_util.time_monitor("cluster_unique_proteins"):
        for row in cur.execute(cluster_query):
            num_read_cluster += 1
            unique_p30s.add(row["p30"])
            unique_p90s.add(row["p90"])
            unique_p100s.add(row["p100"])
        ggdb_logging.info("Done reading clusters")
        ggdb_logging.info(f"num_clusters_rows: {num_read_cluster}")
        ggdb_logging.info(f"num_unique_p30s: {len(unique_p30s)}")
        ggdb_logging.info(f"num_unique_p90s: {len(unique_p90s)}")
        ggdb_logging.info(f"num_unique_p100s: {len(unique_p100s)}")

    # read through neighbor table to gather additional unique p100s
    ggdb_logging.info("Reading prot2protwindow for unique proteins")
    num_prot2protwindow_rows = 0
    with profile_util.memory_monitor("prot2protwindow_unique_proteins"):
        for row in cur.execute(neighbor_query):
            num_prot2protwindow_rows += 1
            unique_p100s.add(row["p1hash"])
            unique_p100s.add(row["p2hash"])

    ggdb_logging.info("Done reading clusters")
    ggdb_logging.info(f"num prot2protwindow_rows: {num_prot2protwindow_rows}")
    ggdb_logging.info(f"num_prot2prot_window_p100s: {2 * num_prot2protwindow_rows}")
    ggdb_logging.info(f"num_unique_p100s: {len(unique_p100s)}")

    ggdb_logging.info("Writing unique proteins to csv files")
    write_neo4j_index_to_file(unique_p30s, p30_fpath, "P30", "p30")
    write_neo4j_index_to_file(unique_p90s, p90_fpath, "P90", "p90")
    write_neo4j_index_to_file(unique_p100s, p100_fpath, "P100", "p100")

    with profile_util.memory_monitor("delete unique sets"):
        del unique_p100s
        del unique_p90s
        del unique_p30s

    unique_p30_p90 = set()
    unique_p90_p100 = set()

    # read through clusters file to gather unique p30_p90, p90_p100
    num_read_cluster_rows = 0
    with profile_util.memory_monitor("clusters_unique_edges"), profile_util.time_monitor("cluster_unique_edges"):
        for row in cur.execute(cluster_query):
            p30_p90 = "_".join([row["p30"], row["p90"]])
            unique_p30_p90.add(p30_p90)
            p90_p100 = "_".join([row["p90"], row["p100"]])
            unique_p90_p100.add(p90_p100)
            num_read_cluster_rows += 1

    ggdb_logging.info(f"num cluster rows: {num_read_cluster_rows}")
    ggdb_logging.info(f"num p30_p90 edges: {len(unique_p30_p90)}")
    ggdb_logging.info(f"num p90_p100 edges: {len(unique_p90_p100)}")

    # write p30_p90 edges, p90_p100 edges to csv
    ggdb_logging.info("Writing cluster_edges to csv files")
    cluster_edges_to_file(unique_p30_p90, "P30", "P90", p30_p90_fpath)
    cluster_edges_to_file(unique_p90_p100, "P90", "P100", p90_p100_fpath)

    # delete edges that we're done with
    with profile_util.memory_monitor("delete p30_90, p90_p100 edges"):
        del unique_p30_p90
        del unique_p90_p100

    # read through neighbor file to gather unique p100_p100s
    with open(p100_neighbor_fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing p100_to_p100_neighbor"):
        csv_fp.write(f":START_ID(P100),:END_ID(P100)\n")
        for row in cur.execute(neighbor_query):
            p1_id = row["p1hash"]
            p2_id = row["p2hash"]
            # sorting unnecessary if we assume these are unique?
            sorted_p1, sorted_p2 = sorted((p1_id, p2_id))
            csv_fp.write(f"{sorted_p1},{sorted_p2}\n")


if __name__ == "__main__":
    main()
