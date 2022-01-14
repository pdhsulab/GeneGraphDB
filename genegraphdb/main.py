import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import sqldb
from genegraphdb import _load
from genegraphdb import _loadsql
from genegraphdb import testing
from genegraphdb import dl_test_data
from genegraphdb import clusternode
from genegraphdb import variables_global as vars_glob
from multiprocessing import Pool, cpu_count
from functools import partial
import os

@click.group()
def cli():
    """A command line tool to create and query a gene graph databases."""
    pass

@cli.command(short_help='Create a new Neo4j GeneGraphDB')
def createdb():
    if not graphdb.hasdb():
        graphdb.createdb()
    else:
        print("Database %s already exists." % DBNAME)

@cli.command(short_help='Get number of nodes and relationships in Neo4j GeneGraphDB')
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

@cli.command(short_help='Clear SQL GeneGraphDB')
def clearsqldb():
    if sqldb.hasdb():
        sqldb.cleardb()
    else:
        print("Database %s does not yet exist." % DBNAME)

@cli.group(short_help='Load data into the graph database.')
def load():
    pass

@load.command(short_help='Load a single sample into the sql database.')
@click.option('--sample-directory', '-s', required=True, help='Path to directory with sample immediately inside.')
@click.option('--sample_id', '-id', required=True, type=str, help='The id of this genome or metagenomic sample.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default="", type=str, help='Any notes on a particular load script runtime')
def singlesql(sample_directory, sample_id, google_bucket, distance, comment):
    if distance is None:
        distance = 5000
    print("start")
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    outfile.close()
    sample_path = sample_directory + sample_id + '/'
    outfilename = sample_path + "ggdb_load_stats.csv"
    _loadsql._single(sample_id, google_bucket, distance, comment, outfilename, sample_directory)
    testing.get_runtime_summarystats(comment, infile_name=sample_path + "/ggdb_load_stats.sql.csv",
                                     outfile_name=sample_path + "/ggdb_summary_stats.sql.csv")

@load.command(short_help='Load multiple samples into the sql database.')
@click.option('--samples-directory', '-s', required=True, help='Path to samples directory. Inside is A001/B001/C001/sampleid/...')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default="", type=str, help='Any notes on a particular load script runtime')
@click.option('--clean_files/--show_temp_files', default=False, help='Remove all temp files generated when loading data into sql database')
def multisql(samples_directory, google_bucket, distance, comment, clean_files):
    if distance is None:
        distance = 5000
    # outfile = open(samples_path + "ggdb_load_stats.csv", "w") #to do - messes with multiprocessing
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    outfile.close()
    #outfilename = samples_path + "ggdb_load_stats.csv"
    outfilename = "ggdb_load_stats.csv"

    # loads samples one by one
    sampleid_to_path_dict = vars_glob.get_sampleid_to_path_dict(samples_directory)
    for key in sampleid_to_path_dict.keys():
        sampleid = key
        samples_path = sampleid_to_path_dict[key]
        if os.path.isdir(samples_path):
            _loadsql._single(sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files)

    testing.get_runtime_summarystats(comment)

    # # Multiprocessing - many samples fail to initially load. these samples are tracked under vars_glob
    # get sample ids and paths to samples. each sample directory contains input data (contigs.tsv, prodigal.faa, prodigal.gff, minced.gff, crisprcastyper.domains.tsv)
    # loadsql_inputs = []
    # for key in vars_glob.drep_samples.keys():
    #     sampleid = key
    #     samples_path = vars_glob.drep_samples[key]
    #     if os.path.isdir(samples_path):
    #         loadsql_inputs.append((sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files))
    # pool = Pool(4)
    # results = pool.starmap(_loadsql._single, loadsql_inputs)
    # pool.close()
    # pool.join()
    
    # # load samples that errored out
    # sample_ids = []
    # for key in vars_glob.drep_samples_error.keys():
    #     sampleid = key
    #     samples_path = vars_glob.drep_samples_error[key]
    #     if os.path.isdir(samples_path):
    #         _loadsql._single(sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files)

#_______________________

@load.command(short_help='Load a single sample into the neo4j database.')
@click.option('--sample_id', '-s', required=True, type=str, help='The id of this genome or metagenomic sample.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default="", type=str, help='Any notes on a particular load script runtime')
def single(sample_id, google_bucket, distance, comment):
    if distance is None:
        distance = 5000
    outfile = open(sample_id + "/ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    _load._single(sample_id, google_bucket, distance, comment, outfile)
    outfile.close()
    testing.get_runtime_summarystats(comment, infile_name=sample_id + "/ggdb_load_stats.csv",
                                     outfile_name=sample_id + "/ggdb_summary_stats.csv")

@load.command(short_help='Load multiple samples into the neo4j database.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default="", type=str, help='Any notes on a particular load script runtime')
@click.option('--clean_files/--show_temp_files', default=False,
              help='Remove all temp files generated when loading data into sql '
                   'database')
def multi(google_bucket, distance, comment, clean_files):
    if distance is None:
        distance = 5000
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    outfile.close()
    outfilename = "ggdb_load_stats.csv"

    
    # get sample ids and paths to samples. each sample directory contains input data (contigs.tsv, prodigal.faa, prodigal.gff, minced.gff, crisprcastyper.domains.tsv)
    load_inputs = []
    for key in vars_glob.drep_samples.keys():
        sampleid = key
        samples_path = vars_glob.drep_samples[key]
        if os.path.isdir(samples_path):
            load_inputs.append(
                (sampleid, google_bucket, distance, comment, outfilename, samples_path, clean_files))
    # # Multiprocessing - has only been successfully tested locally. Neo4j server setup is slow on gcloud.
    pool = Pool(cpu_count())
    results = pool.starmap(_load._single, load_inputs)
    pool.close()
    pool.join()
    # testing.get_runtime_summarystats(comment)
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