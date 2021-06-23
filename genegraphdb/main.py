import click


@click.group()
def cli():
    """A command line tool to create and query a gene graph databases."""
    pass

@cli.command(short_help='Create a new GeneGraphDB')
@click.argument('db')
def createdb(db):
    pass

@cli.group(short_help='Load data into the graph database.')
def load():
    pass

@load.command(short_help='Load a single sample into the database.')
@click.argument('db')
@click.option('--fasta', '-f', required=True, help='The nucleotide FASTA file of the genome.', type=click.Path(exists=True))
@click.option('--protein', '-p', required=True, help='The protein FASTA file.', type=click.Path(exists=True))
@click.option('--gff', '-p', required=True, help='The GFF file.', type=click.Path(exists=True))
@click.option('--google-bucket', '-gb', default=None, help='The Google bucket to store sequences.')
def single(db, fasta, protein, gff, google_bucket):
    print(db, fasta, protein, gff)
    pass

@load.command(short_help='Load multiple samples into the database.')
@click.argument('db')
def multi(db):
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