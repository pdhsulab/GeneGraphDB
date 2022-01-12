from ggdbfetch import clusters, sample2protein, regions
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
    all_p100, p100_to_p90, p100_to_p30 = clusters.get_clusters(target_id, dbpath)
    sample2p100s = sample2protein.get_sample_to_p100s(all_p100, dbpath, sample2path)

    print(sample2p100s)
    for samp in sample2p100s:
        contigs, gene_coords = regions.get_regions(
            sample2path[samp], sample2p100s[samp], p100_to_p90, p100_to_p30, dbpath
        )

        print(contigs)