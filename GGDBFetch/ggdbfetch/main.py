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
def targets_and_baits(infile):
    pass

