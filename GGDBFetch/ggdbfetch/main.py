import click
import os

@click.group()
def cli():
    """A command line tool to fetch ggdb data"""
    pass

@cli.command(short_help='Create a new Neo4j GeneGraphDB')
@click.argument('cluster1')
def createdb():
    if not graphdb.hasdb():
        graphdb.createdb()
    else:
        print("Database %s already exists." % DBNAME)

