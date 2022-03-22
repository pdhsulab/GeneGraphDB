import json
import os

from db_inference.simple_neo4j import SimpleNeo4j
from utils import ggdb_logging
from utils import profile_util

BAITS_PER_FILE = 150
MIN_NUM_CONNECTIONS = 2


def get_colocalization_scores(baits, bait_type):
    neo4j_db = SimpleNeo4j()
    processed_bait_p30s = set()
    p30_to_count_p100s = {}
    colocalization_scores = {}

    # for bait in baits, find neighbors and compute icity for each
    for bait in baits:
        ggdb_logging.info(f"Running bait {bait}")
        bait_p30 = neo4j_db.get_p30_for_protein(bait)
        if bait_p30 in processed_bait_p30s:
            ggdb_logging.info(f"Already computed targets and colocaliztion for bait P30 {bait_p30}")

        if bait_p30 not in p30_to_count_p100s:
            num_p100s = neo4j_db.get_num_p100s(bait_p30)
            p30_to_count_p100s[bait_p30] = num_p100s
        num_bait_p100s = p30_to_count_p100s[bait_p30]
        ggdb_logging.info(f"Bait {bait} has {num_bait_p100s} P100s")

        tgt_p30s_to_conn_counts = neo4j_db.get_targets_and_num_shared_for_bait(bait_p30, MIN_NUM_CONNECTIONS)
        ggdb_logging.info(f"Bait {bait} has {len(tgt_p30s_to_conn_counts)} P30 neighbors with at least  "
                          f"{MIN_NUM_CONNECTIONS} connections")

        # for each target, compute icity (if not cached)
        for tgt_p30, conn_counts in tgt_p30s_to_conn_counts.items():
            tgt_first_key = f"{tgt_p30}|{bait_p30}"
            if tgt_first_key in colocalization_scores:
                ggdb_logging.debug(f"cache hit for {tgt_first_key}")
                continue
            ggdb_logging.debug(f"Computing icity for {tgt_first_key}")

            if tgt_p30 not in p30_to_count_p100s:
                num_p100s = neo4j_db.get_num_p100s(tgt_p30)
                p30_to_count_p100s[tgt_p30] = num_p100s
            num_tgt_p100s = p30_to_count_p100s[tgt_p30]

            if num_tgt_p100s == 0:
                ggdb_logging.debug(f"Skipping nonexistent tgt {tgt_p30}")
                continue

            colocalization_result = {
                "tgt_p30": tgt_p30,
                "bait_p30": bait_p30,
                "num_tgt_p100s": num_tgt_p100s,
                "num_bait_p100s": num_bait_p100s,
                "num_connected_tgt_p100s": conn_counts['num_conn_tgt_p100s'],
                "num_connected_bait_p100s": conn_counts['num_conn_bait_p100s'],
                "num_connections": conn_counts['num_connections'],
                "tgt_colocalization": conn_counts['num_conn_tgt_p100s'] / num_tgt_p100s,
                "bait_colocalization": conn_counts['num_conn_bait_p100s'] / num_bait_p100s,
                "bait_type": bait_type,
            }

            colocalization_scores[tgt_first_key] = colocalization_result
            ggdb_logging.debug(f"added score for {tgt_first_key}")

        processed_bait_p30s.add(bait_p30)

    return colocalization_scores


def main():
    for input_file, bait_type in [
        ("/GeneGraphDB/data/jacob_baits_20220202/cas1.txt", "cas1"),
        ("/GeneGraphDB/data/jacob_baits_20220202/cas2.txt", "cas2"),
        ("/GeneGraphDB/data/jacob_baits_20220202/tnpBs_in_testdb.p100.1e4.txt", "tnpB"),
    ]:
        with profile_util.time_monitor(f"Processing {input_file}"):
            output_file = os.path.join(
                "/GeneGraphDB/data/20220308_neo4j_colocalization/",
                os.path.basename(input_file).replace(".txt", ".json"),
            )

            os.makedirs(os.path.dirname(output_file), exist_ok=True)

            # get baits
            with open(input_file, "r") as f:
                baits = [line.strip() for line in f.readlines()]
            # jacob said cas files accidentally have 20 chars instead of 18
            if "cas" in bait_type:
                baits = [b[:18] for b in baits]

            ggdb_logging.info(f"Found {len(baits)} baits in file {input_file}")

            if BAITS_PER_FILE is not None:
                baits = baits[:BAITS_PER_FILE]
                ggdb_logging.warning(f"Reducing baits to {BAITS_PER_FILE} for debugging purposes")

            colocalization_results = get_colocalization_scores(baits, bait_type)

            # save icity results to a json file for later use
            with open(output_file, "w") as fp:
                json.dump(colocalization_results, fp, indent=2)

            ggdb_logging.info(f"Wrote to file {output_file}")


if __name__ == "__main__":
    main()
