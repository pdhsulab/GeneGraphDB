from ggdbfetch import clusters


def _targets_and_baits(infile, dbpath):

    with open(infile) as inf:

        for line in inf:
            target_id, bait_ids = line.strip().split('\t')
            retrieve_target_and_bait(target_id, bait_ids, dbpath)
            break

def retrieve_target_and_bait(target_id, bait_ids, dbpath):
    all_p100 = clusters.p30_to_p100(target_id, dbpath)
    print(all_p100)