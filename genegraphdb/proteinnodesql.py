import os
import sqlite3
import time
from collections import deque
from csv import reader


def connect_proteins_crisprs(sample_id, max_distance, samples_path = ''):
    sample_id_path = samples_path + sample_id + "/"
    print("Loading protein2protein edges...")
    tic = time.time()
    create_all_protein_crispr_edge_csv(sample_id, max_distance, samples_path)
    load_csv(sample_id_path + "protein2protein.tmp.sql.csv",
             sample_id_path + "protein2crispr.tmp.sql.csv",
             sample_id_path + "protein2protein_window.tmp.sql.csv",
             sample_id_path + "protein2crispr_window.tmp.sql.csv")


    toc = time.time()
    print("Loading protein2protein edges took %f seconds" % (toc - tic))
    return toc-tic

def create_all_protein_crispr_edge_csv(sample_id, max_distance, samples_path = ''):
    sample_id_path = samples_path + sample_id + "/"
    make_merged_coords_csv(sample_id, samples_path)
    outfile_prot_pair = open(sample_id_path + "protein2protein.tmp.sql.csv", "w")
    outfile_prot_crispr_pair = open(sample_id_path + "protein2crispr.tmp.sql.csv", "w")
    outfile_p2p_base_window = open(sample_id_path + "protein2protein_window.tmp.sql.csv", "w")
    outfile_p2c_base_window = open(sample_id_path + "protein2crispr_window.tmp.sql.csv", "w")
    print("recid,phash,qhash", file=outfile_prot_pair)
    print("recid,phash,qhash", file=outfile_prot_crispr_pair)
    print("recid,phash,qhash", file=outfile_p2p_base_window)
    print("recid,phash,qhash", file=outfile_p2c_base_window)

    merge_sorted_coords_csv = sample_id_path + "merged_sorted_coords.tmp.sql.csv" 
        # creates csvs with all pairings between proteins and CRISPRs that are either adjacent ('pairs') or within a genomic neighbourhood ('neighbours')
    create_protein_pair_csv(merge_sorted_coords_csv, outfile_prot_pair, outfile_prot_crispr_pair)
    create_protein_window_csv(merge_sorted_coords_csv, max_distance, outfile_p2p_base_window, outfile_p2c_base_window)

    outfile_prot_pair.close(), outfile_prot_crispr_pair.close(), outfile_p2p_base_window.close(), outfile_p2c_base_window.close()

# allows proteins near CRISPRs to be mapped to each other as neighbours
def make_merged_coords_csv(sample_id, samples_path = ''):
    sample_id_path = samples_path + sample_id + "/"
    os.system("cat " + sample_id_path + "gene_coords.tmp.sql.csv | sed -e '1s/phash/hash/' | cut -d',' -f 1-5 | "
                                        "sed '1s/$/,is_crispr/; 2,$s/$/,0/' > " + sample_id_path + "gene_coords_m.tmp.sql.csv")
    os.system("cat " + sample_id_path + "crispr_coords.tmp.sql.csv | awk 'FNR > 1' | sed '1,$s/$/,1/' > "
              + sample_id_path + "tmp.crispr_coords.sql.csv")
    os.system("cat " + sample_id_path + "gene_coords_m.tmp.sql.csv " + sample_id_path + "tmp.crispr_coords.sql.csv > "
              + sample_id_path + "merged_coords.tmp.sql.csv")
    os.system("head -n1 " + sample_id_path + "merged_coords.tmp.sql.csv > " + sample_id_path +
              "merged_sorted_coords.tmp.sql.csv && tail -n+2 " + sample_id_path + "merged_coords.tmp.sql.csv | sort "
              "--field-separator=',' -k1,1 -k4,4n >> " + sample_id_path + "merged_sorted_coords.tmp.sql.csv")
    os.system("rm " + sample_id_path + "gene_coords_m.tmp.sql.csv " + sample_id_path + "tmp.crispr_coords.sql.csv " + sample_id_path + "merged_coords.tmp.sql.csv")

def create_protein_pair_csv(gene_coords_csv, outfile_prot_pair, outfile_prot_crispr_pair):
    with open(gene_coords_csv, 'r') as g:
        infile = reader(g)
        header = next(infile)
        row1 = next(infile)
        recid, old_phash, old_chash, old_start_coord, old_is_crispr = row1[0], row1[1], row1[2], row1[3], row1[5]
        if header is not None:
            for line in infile:
                recid, cur_phash, cur_chash, cur_start_coord, cur_is_crispr = line[0], line[1], line[2], line[3], line[5]
                if newGene_is_same_contig(old_chash, cur_chash) and is_protein2protein(old_is_crispr, cur_is_crispr):
                    print(recid + "," + cur_phash + "," + old_phash, file=outfile_prot_pair)
                elif newGene_is_same_contig(old_chash, cur_chash):
                    print(recid + "," + cur_phash + "," + old_phash, file=outfile_prot_crispr_pair)
                old_phash = cur_phash
                old_chash = cur_chash
                old_is_crispr = cur_is_crispr

def create_protein_window_csv(merge_coords_csv, max_distance, outfile_p2p, outfile_p2c):
    base_neigh_queue = deque()
    with open(merge_coords_csv, 'r') as f:
        infile = reader(f)
        header = next(infile)
        row1 = next(infile)
        recid, old_phash, old_chash, old_start_coord, old_is_crispr = row1[0], row1[1], row1[2], row1[3], row1[5]
        base_neigh_queue.appendleft({"phash": old_phash, "start_coord": old_start_coord})
        if header is not None:
            for line in infile:

                recid, cur_phash, cur_chash, cur_start_coord, cur_is_crispr = line[0], line[1], line[2], line[3], line[5]
                for _dict in base_neigh_queue:
                    old_phash = _dict["phash"]
                    if newGene_is_same_contig(old_chash, cur_chash) and is_protein2protein(old_is_crispr, cur_is_crispr):
                        print(recid + "," + cur_phash + "," + old_phash, file=outfile_p2p)
                    elif newGene_is_same_contig(old_chash, cur_chash) and is_protein2crispr(old_is_crispr, cur_is_crispr) and int(cur_is_crispr):
                        print(recid + "," + cur_phash + "," + old_phash, file=outfile_p2c)

                base_neigh_queue = update_base_neigh_queue(base_neigh_queue, cur_phash, old_chash,
                                                       cur_chash, max_distance, cur_start_coord,
                                                       old_start_coord)
                old_start_coord = cur_start_coord
                old_chash = cur_chash
                old_is_crispr = cur_is_crispr

# the queue is a sliding window that defines a query protein or crispr's genomic neighbourhood
def update_base_neigh_queue(queue, cur_phash, old_chash, cur_chash, max_distance, new_coord,
                        old_coord):
    if newGene_is_same_contig(old_chash, cur_chash):
        while (int(new_coord) - int(queue[-1]["start_coord"])) > max_distance:
            queue.pop()
            if len(queue) == 0:
                break
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    else:
        queue = deque()
        queue.appendleft({"phash": cur_phash, "start_coord": new_coord})
    return queue

def newGene_is_same_contig(old_chash, cur_chash):
    return old_chash == cur_chash

def is_protein2protein(old_node_is_crispr, cur_node_is_crispr):
    return not (int(old_node_is_crispr) or int(cur_node_is_crispr))

def is_protein2crispr(old_node_is_crispr, cur_node_is_crispr):
    return (int(old_node_is_crispr) and not int(cur_node_is_crispr)) or (not int(old_node_is_crispr) and int(cur_node_is_crispr))

def load_protein_coords(sample_id, samples_path=''):
    con = sqlite3.connect('genegraph.db')
    con.execute("PRAGMA journal_mode=WAL")
    cur = con.cursor()
    gene_coords_csv_path = samples_path + sample_id + '/gene_coords.tmp.sql.csv'
    gene_coords_csv = open(gene_coords_csv_path)
    rows = csv.reader(gene_coords_csv)
    next(rows)
    cmd = '''
        INSERT OR IGNORE INTO proteincoords (phash, contighash, start, end, orientation) VALUES (?,?,?,?,?)
        '''
    cur.executemany(cmd, ((rec[1], rec[2], rec[3], rec[4], rec[5]) for rec in rows))
    con.commit()
    con.close()

# all input csvs have two columns - one for donor protein's phash, other for recipient
def load_csv(csv_path_p2p, csv_path_p2c, csv_path_p2p_window, csv_path_p2c_window):
    con = sqlite3.connect('genegraph.db')
    con.execute("PRAGMA journal_mode=WAL")
    cur = con.cursor()
    p2p_csv, p2c_csv, p2pwindow_csv, p2cwindow_csv= open(csv_path_p2p), open(csv_path_p2c), \
                                                    open(csv_path_p2p_window), open(csv_path_p2c_window)

    p2p_rows = csv.reader(p2p_csv)
    next(p2p_rows)
    p2c_rows = csv.reader(p2c_csv)
    next(p2c_rows)
    p2pwindow_rows = csv.reader(p2pwindow_csv)
    next(p2pwindow_rows)
    p2cwindow_rows = csv.reader(p2cwindow_csv)
    next(p2cwindow_rows)
    p2p_cmd = '''INSERT OR IGNORE INTO prot2prot (p1hash, p2hash) VALUES (?,?)'''
    p2c_cmd = '''INSERT OR IGNORE INTO prot2crispr (p1hash, crisprhash) VALUES (?,?)'''
    p2pwindow_cmd = '''INSERT OR IGNORE INTO prot2protwindow (p1hash, p2hash) VALUES (?,?)'''
    p2cwindow_cmd = '''INSERT OR IGNORE INTO prot2crisprwindow (p1hash, crisprhash) VALUES (?,?)'''

    cur.executemany(p2p_cmd, ((rec[1], rec[2]) for rec in p2p_rows))
    cur.executemany(p2c_cmd, ((rec[1], rec[2]) for rec in p2c_rows))
    cur.executemany(p2pwindow_cmd, ((rec[1], rec[2]) for rec in p2pwindow_rows))
    cur.executemany(p2cwindow_cmd, ((rec[1], rec[2]) for rec in p2cwindow_rows))
    con.commit()
    con.close()


