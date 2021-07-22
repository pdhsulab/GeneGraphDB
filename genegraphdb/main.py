import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load
from genegraphdb import _loadmulti
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
@click.option('--gene-neighbors/--gene-window', default=True, help='Calculate neighbors using number of genes away vs. base pair window.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
def single(sample_id, google_bucket, gene_neighbors, distance, comment):
    # to do - comment out first condition?
    if gene_neighbors and distance is None:
        distance = 3
    elif not gene_neighbors and distance is None:
        distance = 5000
    _load._single(sample_id, google_bucket, gene_neighbors, distance, comment)

@load.command(short_help='Load multiple samples into the database.')
@click.option('--samples_id_path', '-s', required=True, help='The path to directory with genome and metagenomic samples.')
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
@click.option('--gene-neighbors/--gene-window', default=True, help='Calculate neighbors using number of genes away vs. base pair window.')
@click.option('--distance', '-d', default=None, type=int, help='The distance in number of neighbors or base pairs')
@click.option('--comment', '-c', default=None, type=str, help='Any notes on a particular load script runtime')
@click.option('--load_indiv/--load_bulk', default=True, help='Load one or multiple samples with a single csv import')
@click.option('--append_stats/--reset_stats', default=True, help='Reset notes and runtimes of load-script runs')
def multi(samples_id_path, google_bucket, gene_neighbors, distance, comment, load_indiv, append_stats):
    os.chdir(samples_id_path)
    if gene_neighbors and distance is None:
        distance = 3
    elif not gene_neighbors and distance is None:
        distance = 5000

    outfile = open("ggdb_load_stats.csv", "a")
    print("sample_id,load_time,p2p_edge_time,comment", file=outfile)
    if not append_stats:
        outfile.truncate(0)
        print("sample_id,load_time,p2p_edge_time,comment", file=outfile)

    if load_indiv:
        for sample_id in os.listdir():
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _load._single(sample_id, google_bucket, gene_neighbors, distance, comment, outfile)
            except NotADirectoryError:
                print(sample_id + " is not a directory")
    if not load_indiv:
        for sample_id in os.listdir():
            try:
                # to do - implement better way to check if the sample_id is actually a directory
                os.chdir(sample_id)
                os.chdir("..")
                _loadmulti._single(sample_id, google_bucket, gene_neighbors, distance, comment, outfile)
            except NotADirectoryError:
                print(sample_id + " is not a directory")
        outfile.close()
        _loadmulti.bulk_load_protein_crispr_edges(distance, gene_neighbors)
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