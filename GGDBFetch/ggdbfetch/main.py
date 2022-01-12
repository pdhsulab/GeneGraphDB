import click
import os

@click.group()
def cli():
    """A command line tool to fetch ggdb data"""
    pass

@cli.group(short_help='Retrieve stuff.')
def retrieve():
    pass

@retrieve.command(short_help='Retrieve by target and baits.')
@click.argument('infile')
@click.option('--dbpath', default='/home/mdurrant/GeneGraphDB/data/rep_genomes')
def targets_and_baits(infile, dbpath):
    with open(infile) as inf:
        for line in inf:
            target_id, bait_ids = line.strip().split('\t')
            print(target_id, bait_ids)
            break

