import csv
import os
import sqlite3
import unittest
from datetime import datetime

import numpy as np

from genegraphdb import variables_global as vars_glob


# all of these tests were performed locally before implementing the graph database on gcloud
# some tests rely on sample directories that are organized differently from the sample directories on gcloud
class testSQLiteNodeConnections(unittest.TestCase):  # to do - change class name; outdated
    def test_allSQLite_tables_full(self):
        sql_tables = vars_glob.sql_tables
        con = sqlite3.connect("genegraph.db")
        cur = con.cursor()
        for table in sql_tables:
            cmd = "SELECT count(*) FROM {}".format(table)
            cur.execute(cmd)
            rv = cur.fetchall()
            self.assertTrue(rv[0][0] > 0)

    def test_neo4J_SQL_filesizematch(self):
        paths_genomes_annot = []
        for dir in os.listdir("genomes_annot_2/"):
            if dir not in vars_glob.exclude_directories:
                sample_path = "genomes_annot_2/" + dir + "/"
                # to do - use this is_dir check in other files!!
                if os.path.isdir(sample_path):
                    paths_genomes_annot.append(sample_path)
        tempSQLfiles_list = vars_glob.temp_files_list
        allfiles_list = [(file.replace("sql.", ""), file) for file in tempSQLfiles_list]
        for sample_path in paths_genomes_annot:
            for filepair in allfiles_list:
                filepaths = (sample_path + filepair[0], sample_path + filepair[1])
                self.assertEqual(os.path.getsize(filepaths[0]), os.path.getsize(filepaths[1]))

    # mtime describes the time a file was last modified
    def test_all_sqlfiles_updated(self):
        paths_genomes_annot = []
        mtimes_sqlfiles = []
        for dir in os.listdir("genomes_annot/"):
            if dir not in vars_glob.exclude_directories:
                sample_path = "genomes_annot/" + dir + "/"
                # to do - use this is_dir check in other files!!
                if os.path.isdir(sample_path):
                    paths_genomes_annot.append(sample_path)
        tempSQLfiles_list = vars_glob.temp_files_list
        for sample_path in paths_genomes_annot:
            for filename in tempSQLfiles_list:
                filepath = sample_path + filename
                mtimes_sqlfiles.append(os.path.getmtime(filepath))
        mtime_ref = mtimes_sqlfiles[0]
        for mtime in mtimes_sqlfiles:
            self.assertTrue(
                np.abs(mtime - mtime_ref) < 100
            )  # all samples from genomes_annot directory should be loaded in the ggdb in under 100 mtime units

    def test_sqlite_dbsize_stable(self):
        init_db_size = os.path.getsize("genegraph.db")
        os.system("ggdb load multisql -s genomes_annot/ -c '' ")
        final_db_size = os.path.getsize("genegraph.db")
        self.assertEqual(init_db_size, final_db_size)


if __name__ == "__main__":
    unittest.main()


def get_runtime_summarystats(
    comment="", prot_time=0, samples_path="", infile_name="ggdb_load_stats.csv", outfile_name="ggdb_summary_stats.csv"
):
    if samples_path != "":
        infile_name = samples_path + infile_name
        outfile_name = samples_path + outfile_name
    outfile = open(outfile_name, "a")
    with open(infile_name, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        p2p_edge_time, cur_load_time = 0, 0
        for line in f:
            line = line.strip("\n").split(",")
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
                line = line.strip("\n").split(",")
                cur_load_time += float(line[1])
                if p2p_edge_time != "null":
                    p2p_edge_time += float(line[2])
        except StopIteration:
            print("end of csv reached")
        if p2p_edge_time == "null":
            cur_load_time += p2p_edge_time
        print(str(cur_load_time) + "," + str(p2p_edge_time) + "," + comment, file=outfile)


def clean_files(sample_id, samples_path):
    files_to_remove = vars_glob.temp_files_list
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


def log_errors_multi_loadneo4j(samples_path, sample_id, exception):
    outfile = open("ggdb_multineo4j_errorlog.csv", "a")
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print(date_time + "," + samples_path + sample_id + " : " + str(exception), file=outfile)
