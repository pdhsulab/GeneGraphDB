import click
from genegraphdb import *
from genegraphdb import graphdb
from genegraphdb import _load

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
@click.option('--sample_id', '-s', required=True, help='The nucleotide FASTA file of the genome.')
@click.option('--fasta', '-f', required=True, help='The nucleotide FASTA file of the genome.', type=click.Path(exists=True))
@click.option('--protein', '-p', required=True, help='The protein FASTA file.', type=click.Path(exists=True))
@click.option('--gff', '-g', required=True, help='The GFF file.', type=click.Path(exists=True))
@click.option('--contigs', '-g', required=True, help='The contigs file.', type=click.Path(exists=True))
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
def single(sample_id, fasta, protein, gff, contigs, google_bucket):
    _load._single(sample_id, fasta, protein, gff, contigs, google_bucket)

@load.command(short_help='Load multiple samples into the database.')
def multi():
    pass

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
@cli.command(short_help='Find close relatives of query protein seuqences')
@click.argument('db')
@click.argument('protein')
def search(db, protein):
    pass