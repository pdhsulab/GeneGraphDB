import os

from db_inference import simple_sql_db
from utils import ggdb_logging, profile_util

CLUSTER_LIMIT = 100000000  # max is ~200M
NEIGHBOR_LIMIT = 100000000  # max is ~1B


def set_to_neo4j_index(protein_hash_set):
    return {prot_hash: idx + 1 for idx, prot_hash in enumerate(sorted(protein_hash_set))}


def write_neo4j_index_to_file(neo4j_index, fpath, node_name, property_name):
    with open(fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing {fpath}"):
        csv_fp.write(f"id:ID({node_name}),{property_name}\n")
        for protein_hash, idx in sorted(neo4j_index.items(), key=lambda item: item[1]):
            csv_fp.write(f"{str(idx)},{protein_hash}\n")


def cluster_edges_to_file(px_to_py_set, px_neo4j_index, py_neo4j_index, px_node_name, py_node_name, fpath):
    with open(fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing {fpath}"):
        csv_fp.write(f":START_ID({px_node_name}),:END_ID({py_node_name})\n")
        for px_py_hash in px_to_py_set:
            px_hash, py_hash = px_py_hash.split("_")
            px_id = px_neo4j_index[px_hash]
            py_id = py_neo4j_index[py_hash]
            csv_fp.write(f"{px_id},{py_id}\n")


def main():
    ggdb_logging.info(f"CLUSTER_LIMIT is {CLUSTER_LIMIT}")
    ggdb_logging.info(f"NEIGHBOR_LIMIT is {NEIGHBOR_LIMIT}")

    sql_db = simple_sql_db.SimpleSqlDb()
    csv_dir = "/GeneGraphDB/data/neo4j/import"
    p100_fpath = os.path.join(csv_dir, "p100.csv")
    p100_neighbor_fpath = os.path.join(csv_dir, "p100_prot2prot_window_p100.csv")
    p90_fpath = os.path.join(csv_dir, "p90.csv")
    p90_p100_fpath = os.path.join(csv_dir, "p90_cluster_p100.csv")
    p30_fpath = os.path.join(csv_dir, "p30.csv")
    p30_p90_fpath = os.path.join(csv_dir, "p30_cluster_p100.csv")

    unique_p30s = set()
    unique_p90s = set()
    unique_p100s = set()

    # read through clusters file to gather unique p30s, p90s, p100s
    ggdb_logging.info("Reading clusters for unique proteins")
    num_read_cluster = 0
    cur = sql_db.conn.cursor()
    with profile_util.memory_monitor("clusters_unique_proteins"), profile_util.time_monitor("cluster_unique_proteins"):
        for row in cur.execute(f"SELECT * FROM clusters LIMIT {CLUSTER_LIMIT}"):
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
        for row in cur.execute(f"SELECT * FROM prot2protwindow LIMIT {NEIGHBOR_LIMIT}"):
            num_prot2protwindow_rows += 1
            unique_p100s.add(row["p1hash"])
            unique_p100s.add(row["p2hash"])

    ggdb_logging.info("Done reading clusters")
    ggdb_logging.info(f"num prot2protwindow_rows: {num_prot2protwindow_rows}")
    ggdb_logging.info(f"num_prot2prot_window_p100s: {2 * num_prot2protwindow_rows}")
    ggdb_logging.info(f"num_unique_p100s: {len(unique_p100s)}")

    ggdb_logging.info(f"Assigning unique proteins neo4j index")
    with profile_util.memory_monitor("sort unique proteins"), profile_util.time_monitor("sort unique proteins"):
        p100_to_neo4j_idx = set_to_neo4j_index(unique_p100s)
        p90_to_neo4j_idx = set_to_neo4j_index(unique_p90s)
        p30_to_neo4j_idx = set_to_neo4j_index(unique_p30s)

    with profile_util.memory_monitor("delete unique sets"):
        del unique_p100s
        del unique_p90s
        del unique_p30s

    ggdb_logging.info("Writing unique proteins to csv files")
    write_neo4j_index_to_file(p30_to_neo4j_idx, p30_fpath, "P30", "p30")
    write_neo4j_index_to_file(p90_to_neo4j_idx, p90_fpath, "P90", "p90")
    write_neo4j_index_to_file(p100_to_neo4j_idx, p100_fpath, "P100", "p100")

    unique_p30_p90 = set()
    unique_p90_p100 = set()

    # read through clusters file to gather unique p30_p90, p90_p100
    num_read_cluster_rows = 0
    with profile_util.memory_monitor("clusters_unique_edges"), profile_util.time_monitor("cluster_unique_edges"):
        for row in cur.execute(f"SELECT * FROM clusters LIMIT {CLUSTER_LIMIT}"):
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
    cluster_edges_to_file(unique_p30_p90, p30_to_neo4j_idx, p90_to_neo4j_idx, "P30", "P90", p30_p90_fpath)
    cluster_edges_to_file(unique_p90_p100, p90_to_neo4j_idx, p100_to_neo4j_idx, "P90", "P100", p90_p100_fpath)

    # delete edges that we're done with
    with profile_util.memory_monitor("delete p30_90, p90_p100 edges"):
        del unique_p30_p90
        del unique_p90_p100

    # read through neighbor file to gather unique p100_p100s
    with open(p100_neighbor_fpath, "w") as csv_fp, profile_util.time_monitor(f"Writing p100_to_p100_neighbor"):
        csv_fp.write(f":START_ID(P100),:END_ID(P100)\n")
        for row in cur.execute(f"SELECT * FROM prot2protwindow LIMIT {NEIGHBOR_LIMIT}"):
            p1_id = p100_to_neo4j_idx[row["p1hash"]]
            p2_id = p100_to_neo4j_idx[row["p2hash"]]
            # sorting unnecessary if we assume these are unique?
            sorted_p1, sorted_p2 = sorted((p1_id, p2_id))
            csv_fp.write(f"{sorted_p1},{sorted_p2}\n")


if __name__ == "__main__":
    main()
