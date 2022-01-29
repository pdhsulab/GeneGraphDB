import sqlite3
from os.path import basename, join

import pandas as pd
from ggdbfetch import misc
from pyfastx import Fasta


def get_regions(sample_path, p100s, p100_to_p90, p100_to_p30, dbpath):

    sample = basename(sample_path)
    out_fna = join(dbpath, sample_path, sample + ".fna.gz")
    out_coords = join(dbpath, sample_path, sample + ".gene_coords.db")
    out_contigs = join(dbpath, sample_path, sample + ".contigs.db")
    coords = get_gene_coords(sample, out_coords, out_contigs, p100s, p100_to_p90, p100_to_p30)
    contig_range_seqs = dict()
    ranges = coords[
        ["target_p100", "contig", "contig_id", "contig_start", "contig_end", "target_orient"]
    ].drop_duplicates()
    for target_p100, contig, contig_id, contig_start, contig_end, target_orient in zip(
        ranges.target_p100,
        ranges.contig,
        ranges.contig_id,
        ranges.contig_start,
        ranges.contig_end,
        ranges.target_orient,
    ):
        contig_seq = get_contig_seq(out_fna, contig, (contig_start, contig_end), target_orient)
        if "O" in contig_seq or "o" in contig_seq:
            print(sample, contig_id, contig_seq, sep="\t")
        contig_range_seqs[contig_id] = contig_seq

    return contig_range_seqs, coords


def flip_coords(coords, contig_range):
    contig_length = contig_range[1] - contig_range[0]
    out_coords = []
    for p100, start, end, orient in coords:
        start = contig_length - start + 1
        end = contig_length - end + 1
        if orient == "+":
            orient = "-"
        elif orient == "-":
            orient = "+"
        out_coords.append([p100, end, start, orient])
    return out_coords


def get_contig_range(start, end, contig_length):
    cstart, cend = start - 12500, end + 12500
    if cstart < 0:
        cstart = 0
    if cend > contig_length:
        cend = contig_length

    return cstart, cend


def get_contig_seq(fna, contig, contig_range, orient):
    fasta = Fasta(fna)
    seq = str(fasta[contig][contig_range[0] : contig_range[1]]).upper()

    if orient == "-":
        seq = misc.revcomp(seq)

    return seq


def get_gene_coords(sample, gene_coords, contigs_db, p100s, p100_to_p90, p100_to_p30):

    con = sqlite3.connect(gene_coords)
    cur = con.cursor()

    contigs = set()
    p100_coords = set()
    for p100 in p100s:
        cmd = "SELECT p100,contig,start,end,orient FROM coords WHERE p100='{}'".format(p100)
        for p100, contig, start, end, orient in cur.execute(cmd):
            contigs.add(contig)
            p100_coords.add((p100, contig, start, end, orient))
    cur.close()
    con.close()

    contig_lengths = dict()
    con = sqlite3.connect(contigs_db)
    cur = con.cursor()
    for c in contigs:
        cmd = "SELECT length FROM contig_lengths WHERE contig='{}'".format(c)
        contig_lengths[c] = cur.execute(cmd).fetchone()[0]
    cur.close()
    con.close()

    conn = sqlite3.connect(gene_coords)
    cur = conn.cursor()

    all_coords = []
    for p100, contig, start, end, orient in p100_coords:
        contig_range = get_contig_range(start, end, contig_lengths[contig])
        cmd = 'SELECT p100,start,end,orient FROM coords WHERE contig="{c}" AND start >= {s} AND end <= {e}'.format(
            c=contig, s=contig_range[0], e=contig_range[1]
        )
        for this_p100, this_start, this_end, this_strand in cur.execute(cmd):
            this_start -= contig_range[0]
            this_end -= contig_range[0]
            coords = [this_p100, this_start, this_end, this_strand]
            if orient == "-":
                coords = flip_coords([coords], contig_range)[0]
            contig_id = (
                p100_to_p30[p100]
                + "|"
                + p100_to_p90[p100]
                + "|"
                + p100
                + "|"
                + sample
                + "|"
                + contig
                + ":"
                + str(contig_range[0])
                + "-"
                + str(contig_range[1])
            )
            if orient == "+":
                contig_id += "(FWD)"
            else:
                contig_id += "(REV)"
            all_coords.append(
                coords + [p100, orient, contig, contig_id, contig_range[0], contig_range[1], contig_lengths[contig]]
            )
    cur.close()
    conn.close()

    coords = pd.DataFrame(all_coords)
    coords.columns = [
        "p100",
        "start",
        "end",
        "strand",
        "target_p100",
        "target_orient",
        "contig",
        "contig_id",
        "contig_start",
        "contig_end",
        "contig_length",
    ]

    return coords
