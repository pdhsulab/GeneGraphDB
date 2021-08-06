import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from genegraphdb import _loadmulti
from genegraphdb import testing
from genegraphdb import dl_test_data
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


@cli.group(short_help='Load data into the graph database.')
def load():
    pass

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
    try:
        os.chdir(samples_id_path)
    except:
        # all test data is one directory up
        os.chdir("..")
        test_data_dir = samples_id_path.replace("../", "")
        dl_test_data.download_dir(test_data_dir)
        os.chdir(test_data_dir)
    if distance is None:
        distance = 5000
    outfile = open("ggdb_load_stats.csv", "w")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    if load_indiv:
        for sample_id in os.listdir():
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _load._single(sample_id, google_bucket, distance, comment, outfile)
            except NotADirectoryError:
                print(sample_id + " is not a directory")
        outfile.close()
        testing.get_runtime_summarystats(comment)
    elif not load_indiv:
        for sample_id in os.listdir():
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _loadmulti._single(sample_id, google_bucket, distance, comment, outfile)
            except NotADirectoryError:
                print(sample_id + " is not a directory")
        outfile.close()
        prot_edge_load_time = _loadmulti.bulk_connect_proteins_crisprs(distance)
        testing.get_runtime_summarystats(comment, prot_edge_load_time)
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