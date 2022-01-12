from ggdbfetch import clusters
from ggdbfetch import sample2protein


def _targets_and_baits(infile, dbpath):

    with open(infile) as inf:

        for line in inf:
            target_id, bait_ids = line.strip().split('\t')
            retrieve_target_and_bait(target_id, bait_ids, dbpath)
            break

def retrieve_target_and_bait(target_id, bait_ids, dbpath):
    all_p100 = clusters.p30_to_p100(target_id, dbpath)
    sample2p100s = sample2protein.get_sample_to_p100s(all_p100, dbpath)

    print(sample2p100s)