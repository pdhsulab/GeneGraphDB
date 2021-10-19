from genegraphdb import *
import csv
import os
from datetime import datetime
import unittest
import pandas as pd
from genegraphdb import variables_global as vars_glob
import sqlite3
import numpy as np

class testSQLiteNodeConnections(unittest.TestCase): # to do - change class name; outdated
    #All tables are populated
    # Neo4J file sizes match SQL file sizes
    # All files are updated properly
    # (After indexing)
        # Number of rows per table is stable between load runs
        # Something duplicate related
    def test_allSQLite_tables_full(self):
        sql_tables = vars_glob.sql_tables
        con = sqlite3.connect('genegraph.db')
        cur = con.cursor()
        for table in sql_tables:
            cmd = "SELECT count(*) FROM {}".format(table)
            cur.execute(cmd)
            rv = cur.fetchall()
            # print(rv,rv[0][0])
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
        allfiles_list = [(file.replace("sql.",""), file) for file in tempSQLfiles_list]
        for sample_path in paths_genomes_annot:
            for filepair in allfiles_list:
                filepaths = (sample_path + filepair[0], sample_path + filepair[1])
                self.assertEqual(os.path.getsize(filepaths[0]), os.path.getsize(filepaths[1]))
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
            # print(mtime, mtime_ref, mtime-mtime_ref) # to do - delete
            self.assertTrue(np.abs(mtime - mtime_ref) < 100) # ~100 mtime units is the around the observed time it takes to load samples in SQLite
    def test_sqlite_dbsize_stable(self):
        init_db_size = os.path.getsize("genegraph.db")
        os.system("ggdb load multisql -s genomes_annot/ -c '' ")
        final_db_size = os.path.getsize("genegraph.db")
        self.assertEqual(init_db_size, final_db_size)

    # def test_num_prot2prot(self):
    #     num_prot2prot_csv =
    #
    #     num_prot2prot_sql =
    def test_example_queries(self):
        querysuccess = True
        con = sqlite3.connect('genegraph.db')
        cur = con.cursor()
        # *optional* find contigs with many proteins (single protein contigs won't join with prot2prot)
        # find a protein's neighbours based on protein's hashid
        # account for duplicate proteins in other contigs
        query_prot2prot = """
        SELECT * FROM proteins as p
        INNER JOIN prot2prot as p2p
        ON (p.hashid = p2p.p1hash OR p.hashid = p2p.p2hash)
        WHERE p.hashid = '188f3c7cbdf9c81f121f|'
        """
        query_prot2crispr = """
        """
        try:
            cur.execute(query_prot2prot)
            #cur.execute(query_prot2crispr)
        except:
            querysuccess = False
        #rv = cur.fetchall()
        con.close()
        self.assertTrue(querysuccess)
        pass

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