import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import sqldb
from genegraphdb import _load
from genegraphdb import _loadsql
from genegraphdb import _loadmulti
# from genegraphdb import _loadmultisql
from genegraphdb import testing
from genegraphdb import dl_test_data
from genegraphdb import clusternode
from genegraphdb import variables_global as vars_glob
from genegraphdb import indexsql
from multiprocessing import Pool, cpu_count
from functools import partial
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
# @click.option('--samples_path', '-s', required=True, help='The path to directory with genome and metagenomic samples. '
#                                                           'Must include / at the end')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
@click.option('--clean_files/--show_temp_files', default=False, help='Remove all temp files generated when loading data into sql database')
def multisql(google_bucket, distance, comment, clean_files):
    if distance is None:
        distance = 5000
    # outfile = open(samples_path + "ggdb_load_stats.csv", "w") #to do - messes with multiprocessing
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    outfile.close()
    #outfilename = samples_path + "ggdb_load_stats.csv"
    outfilename = "ggdb_load_stats.csv"
    sample_ids = []
    
    for key in vars_glob.drep_samples_error.keys():
        sampleid = key
        samples_path = vars_glob.drep_samples_error[key]
        if os.path.isdir(samples_path):
            _loadsql._single(sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files)

    # # Code that works for different sample file directory structure - genomes_annot
    # for sample_id in os.listdir(samples_path):
    #     if os.path.isdir(samples_path + sample_id):
    #         sample_ids.append(sample_id)
    #         # try:
    #         #     _loadsql._single(sample_id, google_bucket, distance, comment, outfile, samples_path, clean_files)
    #         # except Exception as e:
    #         #     testing.log_errors_multisql_loadsql(samples_path, sample_id, e)
    # loadsql_inputs = [(sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files)
    #                   for sampleid in sample_ids]

    # #this for loop reloads all files that failed to previously load
    # for key in vars_glob.drep_samples_error.keys():
    #     sampleid = key
    #     samples_path = vars_glob.drep_samples_error[key]
    #     if os.path.isdir(samples_path):
    #         _loadsql._single(sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files)
    
    # Multiprocessing - leads to many errors; to do - figure out how to resolve this
    # loadsql_inputs = []
    # print(vars_glob.drep_samples)
    # for key in vars_glob.drep_samples.keys():
    #     sampleid = key
    #     samples_path = vars_glob.drep_samples[key]
    #     if os.path.isdir(samples_path):
    #         loadsql_inputs.append((sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files))
    # print(loadsql_inputs)
    # pool = Pool(cpu_count())
    # results = pool.starmap(_loadsql._single, loadsql_inputs)
    # pool.close()
    # pool.join()

    #outfile.close()
    #testing.get_runtime_summarystats(comment, samples_path=samples_path)
    testing.get_runtime_summarystats(comment)

    # clusternode.load_cluster_nodes() #to do - SQLyfy
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
# @click.option('--samples_path', '-s', required=True, help='The path to directory with genome and metagenomic samples. Must end with '
#                                                              '/')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
@click.option('--clean_files/--show_temp_files', default=False,
              help='Remove all temp files generated when loading data into sql '
                   'database')
def multi(google_bucket, distance, comment, clean_files):
    if distance is None:
        distance = 5000
    # outfile = open(samples_path + "ggdb_load_stats.csv", "w") # to do - messes w/ multiprocessing
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    outfile.close()
    outfilename = "ggdb_load_stats.csv"

    # Multiprocessing - to do - see if this results in weird errors
    load_inputs = []
    print(vars_glob.drep_samples)
    for key in vars_glob.drep_samples.keys():
        sampleid = key
        samples_path = vars_glob.drep_samples[key]
        if os.path.isdir(samples_path):
            load_inputs.append(
                (sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files))
    print(load_inputs)
    pool = Pool(cpu_count())
    results = pool.starmap(_load._single, load_inputs)
    pool.close()
    pool.join()
    #outfile.close()
    #testing.get_runtime_summarystats(comment, samples_path=samples_path)

    # clusternode.load_cluster_nodes()

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