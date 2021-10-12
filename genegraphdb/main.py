import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import sqldb
from genegraphdb import _load
from genegraphdb import _loadsql
from genegraphdb import _loadmulti
from genegraphdb import _loadmultisql
from genegraphdb import testing
from genegraphdb import dl_test_data
from genegraphdb import clusternode
from genegraphdb import variables_global as vars_glob
import os

@click.group()
def cli():
    """A command line tool to create and query a gene graph databases."""
    pass

@cli.command(short_help='Create a new GeneGraphDB')
def createdb():

    if not graphdb.hasdb():
        graphdb.createdb()
    else:
        print("Database %s already exists." % DBNAME)

@cli.command(short_help='Create a new GeneGraphDB')
def dbinfo():

    if not graphdb.hasdb():
        graphdb.createdb()

    print("Total nodes in database:", graphdb.num_nodes())
    print("Total relationships in database:", graphdb.num_rels())

@cli.command(short_help='Create a new SQL GeneGraphDB')
def createsqldb():
    if not sqldb.hasdb():
        sqldb.createdb()
    else:
        print("Database %s already exists." % DBNAME)

@cli.command(short_help='Create a new SQL GeneGraphDB')
def sqldbinfo():
    if not sqldb.hasdb():
        sqldb.createdb()
    print("Total nodes in database:", graphdb.num_nodes())
    print("Total relationships in database:", graphdb.num_rels())

@cli.command(short_help='Create a new SQL GeneGraphDB')
def clearsqldb():
    if sqldb.hasdb():
        sqldb.cleardb()
    else:
        print("Database %s does not yet exist." % DBNAME)

@cli.group(short_help='Load data into the graph database.')
def load():
    pass

#_______________________
@load.command(short_help='Load a single sample into the sql database.')
@click.option('--sample_id', '-s', required=True, type=str, help='The id of this genome or metagenomic sample.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
def singlesql(sample_id, google_bucket, distance, comment):
    # to do - comment out first condition?
    if distance is None:
        distance = 5000
    outfile = open(sample_id + "/ggdb_load_stats.sql.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    _loadsql._single(sample_id, google_bucket, distance, comment, outfile)
    outfile.close()
    testing.get_runtime_summarystats(comment, infile_name=sample_id + "/ggdb_load_stats.sql.csv",
                                     outfile_name=sample_id + "/ggdb_summary_stats.sql.csv")

@load.command(short_help='Load multiple samples into the sqldatabase.')
@click.option('--samples_path', '-s', required=True, help='The path to directory with genome and metagenomic samples. '
                                                          'Must include / at the end')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
@click.option('--load_indiv/--load_bulk', default=True, help='Load one or multiple samples with a single csv import')
@click.option('--clean_files/--show_temp_files', default=False, help='Remove all temp files generated when loading data into sql '
                                                         'database')
def multisql(samples_path, google_bucket, distance, comment, load_indiv, clean_files):
    # try:
    #     os.chdir(samples_path)
    # except:
    #     # all test data is one directory up
    #     os.chdir("..")
    #     test_data_dir = samples_path.replace("../", "")
    #     dl_test_data.download_dir(test_data_dir)
    #     os.chdir(test_data_dir)
    if distance is None:
        distance = 5000
    outfile = open(samples_path + "ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    if load_indiv:
        for sample_id in os.listdir(samples_path):
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(samples_path + sample_id)
                os.chdir("../..")
                try:
                    _loadsql._single(sample_id, google_bucket, distance, comment, outfile, samples_path, clean_files)
                except Exception as e:
                    testing.log_errors_multisql_loadsql(samples_path, sample_id, e)
            except NotADirectoryError:
                pass
        outfile.close()
        testing.get_runtime_summarystats(comment, samples_path=samples_path)
    elif not load_indiv:
        for sample_id in os.listdir(samples_path):
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(samples_path + sample_id)
                os.chdir("../..")
                _loadmultisql._single(sample_id, google_bucket, distance, comment, outfile, samples_path, clean_files)
            except NotADirectoryError:
                pass
        outfile.close()
        prot_edge_load_time = _loadmultisql.bulk_connect_proteins_crisprs(distance, samples_path)
        testing.get_runtime_summarystats(comment, prot_time=prot_edge_load_time, samples_path=samples_path)
    os.chdir("..")
#_______________________

@load.command(short_help='Load a single sample into the database.')
@click.option('--sample_id', '-s', required=True, type=str, help='The id of this genome or metagenomic sample.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
def single(sample_id, google_bucket, distance, comment):
    # to do - comment out first condition?
    if distance is None:
        distance = 5000
    outfile = open(sample_id + "/ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    _load._single(sample_id, google_bucket, distance, comment, outfile)
    outfile.close()
    testing.get_runtime_summarystats(comment, infile_name=sample_id + "/ggdb_load_stats.csv",
                                     outfile_name=sample_id + "/ggdb_summary_stats.csv")

@load.command(short_help='Load multiple samples into the database.')
@click.option('--samples_id_path', '-s', required=True, help='The path to directory with genome and metagenomic samples.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
@click.option('--load_indiv/--load_bulk', default=True, help='Load one or multiple samples with a single csv import')
def multi(samples_id_path, google_bucket, distance, comment, load_indiv):
    # try:
    #     os.chdir(samples_id_path)
    # except:
    #     # all test data is one directory up
    #     os.chdir("..")
    #     test_data_dir = samples_id_path.replace("../", "")
    #     dl_test_data.download_dir(test_data_dir)
    #     os.chdir(test_data_dir)
    if distance is None:
        distance = 5000
    outfile = open(samples_id_path + "/ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    if load_indiv:
        for sample_id in os.listdir():
            if sample_id in vars_glob.exclude_directories:
                print(sample_id)
                continue
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _load._single(sample_id, google_bucket, distance, comment, outfile)
            except NotADirectoryError:
                pass
        outfile.close()
        testing.get_runtime_summarystats(comment, samples_path=samples_path)
    elif not load_indiv:
        for sample_id in os.listdir():
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _loadmulti._single(sample_id, google_bucket, distance, comment, outfile)
            except NotADirectoryError:
                pass
        outfile.close()
        prot_edge_load_time = _loadmulti.bulk_connect_proteins_crisprs(distance)
        testing.get_runtime_summarystats(comment, prot_time=prot_edge_load_time, samples_path=samples_path)
    clusternode.load_cluster_nodes()
    os.chdir("..")

@cli.command(short_help='Ådd clusters to protein nodes in the database.')
def addclusters():
    pass

@cli.command(short_help='Get the association score for two gene clusters.')
@click.argument('db')
@click.argument('cluster1')
@click.argument('cluster2')
def association(db, cluster1, cluster2):
    pass

# This is Sam's job
@cli.command(short_help = 'Prepare database for queries')
def prepsearchdb():
    graphdb.kmerdb()
    pass

@cli.command(short_help='Find close relatives of query protein seuqences')
@click.argument('db')
@click.argument('protein')
def search(db, protein):
    graphdb.search(protein)
    pass