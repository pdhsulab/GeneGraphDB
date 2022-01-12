from ggdbfetch import clusters
from ggdbfetch import sample2protein
from ggdbfetch import SAMPLE2PATH
from os.path import join

def _targets_and_baits(infile, dbpath):

    print("Loading sample2path")
    sample2path = dict()
    with open(join(dbpath, SAMPLE2PATH)) as inf:
        for line in inf:
            sample, spath = line.strip().split()
            sample2path[sample] = spath


    with open(infile) as inf:
        for line in inf:
            target_id, bait_ids = line.strip().split('\t')
            retrieve_target_and_bait(target_id, bait_ids, dbpath, sample2path)
            break

def retrieve_target_and_bait(target_id, bait_ids, dbpath, sample2path):
    all_p100 = clusters.p30_to_p100(target_id, dbpath)
    sample2p100s = sample2protein.get_sample_to_p100s(all_p100, dbpath, sample2path)

    for samp in sample2p100s:
        print(samp, len(sample2p100s[samp]))