from collections import defaultdict
from typing import Dict, Tuple

import networkx as nx

from db_inference.simple_sql_db import SimpleSqlDb

# hex colors (for better graphing in pyvis)
BLACK = "#000000"
WHITE = "#FFFFFF"
RED = "#FF0000"
GREEN = "#00FF00"
BLUE = "#0000FF"
YELLOW = "#FFFF00"
CYAN = "#00FFFF"
MAGENTA = "#FF00FF"
PURPLE = "#A020F0"


def build_icity_graph(sql_db: SimpleSqlDb, tgt_p30: str, bait_p30: str) -> nx.Graph:
    """Builds a networkx graph that makes calculating icity relatively clean.  Simplifies visualizing protein."""
    G = nx.Graph()
    p30_to_p100s = defaultdict(set)
    # add target cluster, bait cluster
    for p30_hash, p100_color, label in ((tgt_p30, BLUE, "target"), (bait_p30, PURPLE, "bait")):
        # add central p30 node
        p30 = p30_hash + "_p30"
        G.add_node(p30, size=50, color=RED)
        for cluster_row in sql_db.get_p30_cluster_members(p30_hash):
            # add other p30 -> p90 connection
            p90 = cluster_row["p90"] + "_p90"
            if not G.has_node(p90):
                G.add_node(p90, size=25, color=YELLOW)
                G.add_edge(p30, p90, type="p30_clustering")
            # add p90 -> p100 connection
            p100 = cluster_row["p100"]
            G.add_node(p100, size=10, color=p100_color)
            G.add_edge(p90, p100, type="p90_clustering")
            p30_to_p100s[p30_hash].add(p100)

    # add windowed neighbor edges
    # NOTE(john): we run a sql query for every p100 reachable from either bait or target.
    # choose the smaller of these two to reduce the number of SQL queries
    query_p30 = min(p30_to_p100s, key=lambda k: len(p30_to_p100s[k]))
    non_query_p30 = bait_p30 if query_p30 == tgt_p30 else tgt_p30
    for query_p100 in p30_to_p100s[query_p30]:
        p100_neighbors = sql_db.get_p100_windowed_neighbors(query_p100)
        # only add an edge if the neighbor belongs to opposite p30's cluster (vs. a p30 that is neither target nor bait)
        for neighbor in p100_neighbors:
            if neighbor in p30_to_p100s[non_query_p30]:
                G.add_edge(query_p100, neighbor)

    # TODO(john): try collapsing above into a single SQL query!
    # profile difference in runtime, while asserting calculated icity values stay the same for SQL heavy implementation

    return G


def compute_icity_on_graph(G: nx.Graph, tgt_p30_hash: str) -> Dict:
    icity_res = {"tgt_hash": tgt_p30_hash}

    # how many target p100's?
    tgt_p30_node = tgt_p30_hash + "_p30"
    tgt_p100s = list(nx.descendants_at_distance(G, tgt_p30_node, 2))
    icity_res["num_tgt_p100s"] = len(tgt_p100s)

    # how many target p100's connected to bait?
    p100_to_bool_icity = {}
    for p100 in tgt_p100s:
        # having a neighbor besides your p90 clustering implies you're connected to a bait p100
        p100_to_bool_icity[p100] = any([not n.endswith("_p90") for n in G.neighbors(p100)])
    num_p100_positive = sum(p100_to_bool_icity.values())
    icity_res["num_tgt_p100s_connected"] = num_p100_positive
    icity_res["icity_p100"] = num_p100_positive / len(tgt_p100s)

    # how many p90's?
    tgt_p90s = list(G.neighbors(tgt_p30_node))
    icity_res["num_tgt_p90s"] = len(tgt_p90s)

    # calculate target p90's w/ connection two different ways
    p90s_any_icity = 0
    p90s_majority_icity = 0
    for p90_cluster in tgt_p90s:
        p100_neighbors = [n for n in G.neighbors(p90_cluster) if n != tgt_p30_node]
        p100_positive = [n for n in p100_neighbors if p100_to_bool_icity[n]]
        if len(p100_positive) > 0:
            p90s_any_icity += 1
        if (len(p100_positive) / len(p100_neighbors)) > 0.5:
            p90s_majority_icity += 1

    icity_res["num_tgt_p90s"] = len(tgt_p90s)
    icity_res["num_tgt_p90s_any_p100_connected"] = p90s_any_icity
    icity_res["num_tgt_p90s_majority_p100_connected"] = p90s_majority_icity

    icity_res["icity_p90_any"] = p90s_any_icity / len(tgt_p90s)
    icity_res["icity p90 majority"] = p90s_majority_icity / len(tgt_p90s)
    return icity_res


def calc_icity_tgt_bait(sql_db: SimpleSqlDb, tgt_p30: str, bait_p30: str) -> Tuple[Dict, Dict, nx.Graph]:
    G = build_icity_graph(sql_db, tgt_p30, bait_p30)
    tgt_res = compute_icity_on_graph(G, tgt_p30)
    bait_res = compute_icity_on_graph(G, bait_p30)
    return tgt_res, bait_res, G
