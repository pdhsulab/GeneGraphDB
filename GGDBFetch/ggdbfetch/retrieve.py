from ggdbfetch import clusters, sample2protein, regions
from ggdbfetch import SAMPLE2PATH
from os.path import join
from multiprocessing import Pool

def _targets_and_baits(infile, dbpath, outdir, threads):

    print("Loading sample2path")
    sample2path = dict()
    with open(join(dbpath, SAMPLE2PATH)) as inf:
        for line in inf:
            sample, spath = line.strip().split()
            sample2path[sample] = spath


    with open(infile) as inf:
        for line in inf:
            target_id, bait_ids = line.strip().split('\t')
            print("Working on", target_id)
            retrieve_target_and_bait(target_id, bait_ids, dbpath, sample2path, outdir, threads)

def retrieve_target_and_bait(target_id, bait_ids, dbpath, sample2path, outdir, threads):

    print('\tGetting clusters...')
    all_p100, p100_to_p90, p100_to_p30 = clusters.get_clusters(target_id, dbpath)
    print('\tGetting samples...')
    sample2p100s = sample2protein.get_sample_to_p100s(all_p100, dbpath, sample2path)

    print('\tGetting regions...')
    args = [(sample2path[samp], sample2p100s[samp], p100_to_p90, p100_to_p30, dbpath) for samp in sample2p100s]
    with Pool(threads) as pool:
        results = pool.starmap(regions.get_regions, args)

    out_gff = open(join(outdir, target_id + '.examples.gff'), 'w')
    for contigs, coords in results:
        for p100, start, end, strand, target_p100, contig in zip(
                coords.p100, coords.start, coords.end, coords.strand, coords.target_p100, coords.contig_id
        ):
            seqname, source, feature, score, frame = contig, "PRODIGAL", "CDS", 1, 0
            if p100 == target_p100:
                feature = 'CDS_target'
            elif p100 in bait_ids:
                feature = 'CDS_bait'
            attrib = 'ID=' + p100
            out = [seqname, source, feature, start, end, score, strand, frame, attrib]
            print(*out, sep='\t', file=out_gff)
    for contigs, coords in results:
        for contig in contigs:
            print('>' + contig, contigs[contig], sep='\n', file=out_gff)
    out_gff.close()