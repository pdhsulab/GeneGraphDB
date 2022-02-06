import json
import os

from db_inference import calc_icity
from db_inference.simple_sql_db import SimpleSqlDb
from utils import ggdb_logging

INPUT_FILE = "/GeneGraphDB/data/jacob_baits_20220202/cas1.txt"
OUTPUT_FILE = os.path.join("/GeneGraphDB/data/icity_results/",
                           os.path.basename(INPUT_FILE).replace(".txt", ".json"))


def main():
    sql_db = SimpleSqlDb()
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    # get baits
    with open(INPUT_FILE, "r") as f:
        baits = [line.strip() for line in f.readlines()]
    # jacob said files here accidentally have 20 chars instead of 18
    baits = [b[:18] for b in baits]
    ggdb_logging.info(f"Found {len(baits)} baits in file {INPUT_FILE}")

    # REMOVE REMOVE (get rid of this limit)
    baits = baits[:10]

    p100_to_p30 = {}  # cache of cluster lookups
    icity_results = {}  # cache of icity results by p30 cluster

    # for bait in baits, find neighbors and compute icity for each
    for bait in baits:
        ggdb_logging.info(f"Running bait {bait}")
        if bait not in p100_to_p30:
            bait_p30 = sql_db.get_p30_cluster_for_p100(bait)["p30"]
            p100_to_p30[bait] = bait_p30
        bait_p30 = p100_to_p30[bait]

        bait_neighbors = sql_db.get_p100_windowed_neighbors(bait)
        ggdb_logging.info(f"Bait {bait} has {len(bait_neighbors)} neighbors")

        # for each target, compute icity (if not cached)
        for tgt in bait_neighbors:
            if tgt not in p100_to_p30:
                tgt_p30_row = sql_db.get_p30_cluster_for_p100(tgt)
                if tgt_p30_row is None:
                    ggdb_logging.info("Skipping missing target p30")
                    continue
                p100_to_p30[tgt] = tgt_p30_row["p30"]
            tgt_p30 = p100_to_p30[tgt]

            tgt_first_key = f"{tgt_p30}|{bait_p30}"
            if tgt_first_key in icity_results:
                ggdb_logging.info(f"cache hit for {tgt_first_key}")
                continue
            ggdb_logging.info(f"Computing icity for {tgt_first_key}")

            icity_graph = calc_icity.build_icity_graph(sql_db, tgt_p30, bait_p30)
            # compute icity in both directions; i expect this will be informative
            tgt_first_icity = calc_icity.compute_icity_on_graph(icity_graph, tgt_p30)
            tgt_first_icity["bait_hash"] = bait_p30
            icity_results[tgt_first_key] = tgt_first_icity

            bait_first_key = f"{bait_p30}|{tgt_p30}"
            bait_first_icity = calc_icity.compute_icity_on_graph(icity_graph, bait_p30)
            bait_first_icity["bait_hash"] = tgt_p30
            icity_results[bait_first_key] = bait_first_icity

    # save icity results to a json file for later use
    with open(OUTPUT_FILE, 'w',) as fp:
        json.dump(icity_results, fp, indent=2)

    ggdb_logging.info(f"Wrote to file {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
