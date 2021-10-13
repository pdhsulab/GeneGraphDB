from genegraphdb import *
import csv
import os
from datetime import datetime
import unittest
import pandas as pd

class testSQLiteNodeConnections(unittest.TestCase):
    def test_100crispr2prot(self):
        pass
    def test_100prot2prot(self):
        # number edges per protein = 2
        # protein window in protein2protein_window.tmp.sql.csv matches
        test_prot2prot_path = "genomes_annot/561291/protein2protein.tmp.sql.csv"
        prot2prot_df = pd.read_csv(test_prot2prot_path).head(500)
        phash_s = prot2prot_df["phash"]
        print(type(phash_s))
        self.assertTrue(phash_s.is_unique)

    def test_protein2protein(self):
        self.assertEqual(1, 1)
        self.assertEqual(1, 1)
        self.assertEqual(1, 1)
if __name__ == '__main__':
    unittest.main()

# To do: ordering of arguments is super confusing: fix
def get_runtime_summarystats(comment="", prot_time = 0, samples_path = "", infile_name = "ggdb_load_stats.csv", outfile_name = "ggdb_summary_stats.csv"):
    if samples_path != "":
        infile_name = samples_path + infile_name
        outfile_name = samples_path + outfile_name
    outfile = open(outfile_name, "a")
    with open(infile_name, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        p2p_edge_time, cur_load_time = 0, 0
        for line in f:
            line = line.strip('\n').split(',')
            cur_load_time = float(line[1])
            p2p_edge_time = line[2]
            if p2p_edge_time != "null":
                p2p_edge_time = float(p2p_edge_time)
            else:
                p2p_edge_time = prot_time
            break
        try:
            next(reader)
            for line in f:
                line = line.strip('\n').split(',')
                cur_load_time += float(line[1])
                if p2p_edge_time != "null":
                    p2p_edge_time += float(line[2])
        except StopIteration:
            print("end of csv reached")
        if p2p_edge_time == "null":
            cur_load_time += p2p_edge_time
        print(str(cur_load_time) + "," + str(p2p_edge_time) + "," + comment, file=outfile)

def clean_files(sample_id, samples_path):
    files_to_remove = ["contig2sample.tmp.sql.csv", "crispr_coords.tmp.sql.csv", "CRISPRs.tmp.sql.csv",
                       "gene_coords.tmp.sql.csv", "merged_sorted_coords.tmp.sql.csv", "merged.sorted.tmp.sql.gff",
                       "protein2crispr_window.tmp.sql.csv", "protein2crispr.tmp.sql.csv", "protein2protein_window.tmp.sql.csv",
                       "protein2protein.tmp.sql.csv", "proteins.tmp.sql.csv", "temp.minced.sql.gff"]
    cmd = "rm "
    for file in files_to_remove:
        cmd += samples_path + sample_id + "/" + file + " "
    print(cmd)
    os.system(cmd)

def log_errors_multisql_loadsql(samples_path, sample_id, exception):
    outfile = open("ggdb_multisql_errorlog.csv", "a")
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print(date_time + "," + samples_path + sample_id + " : " + str(exception), file=outfile)