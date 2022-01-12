import click
import os
from ggdbfetch.retrieve import _targets_and_baits
import shutil
import sys

@click.group()
def cli():
    """A command line tool to fetch ggdb data"""
    pass

@cli.group(short_help='Retrieve stuff.')
def retrieve():
    pass

@retrieve.command(short_help='Retrieve by target and baits.')
@click.argument('infile')
@click.option('--outdir', '-o', default='ggdbfetch_output')
@click.option('--dbpath', default='/home/mdurrant/GeneGraphDB/data/rep_genomes')
@click.option('--force/--no-force', default=False)
@click.option('--threads', '-t', default=1)
def targets_and_baits(infile, outdir, dbpath, force, threads):

    if os.path.isdir(outdir) and not force:
        print("Output directory already exists, exiting...")
        sys.exit()
    elif os.path.isdir(outdir):
        print("Output directory exists, but removing because --force flag is given...")
        shutil.rmtree(outdir)

    os.makedirs(outdir)
    _targets_and_baits(infile, dbpath, outdir, threads)


