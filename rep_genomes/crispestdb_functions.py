from google.cloud import storage
import numpy as np
import pyfastx
import sqlite3
import subprocess
from os import remove, devnull
from itertools import cycle
from multiprocessing import Pool
import math
from pandas.errors import EmptyDataError
from collections import defaultdict
import timeit
import csv
from snakemake import shell
from snakemake.io import InputFiles, OutputFiles, Params
from pathlib import Path
import shutil
import re
import gzip
from Bio import SeqIO
import os
from os.path import basename, join, dirname, isfile
from os import rename, makedirs
import wget
import psycopg2
import hashlib
import sys
from glob import glob
from random import randint
import pandas as pd

psql_user, psql_pass, psql_host, psql_db = 'postgres', '3c28MP19clMfd5HP', '10.27.128.7', 'crispestdb'
pwd = os.path.dirname(os.path.realpath(__file__))
home = str(Path.home())

def remove(f):
	if os.path.exists(f):
		os.remove(f)
	else:
		pass

def get_connection():
	conn = psycopg2.connect(user=psql_user, password=psql_pass, host=psql_host, database=psql_db)
	return conn

def get_sample_field(field, params, conn):
	with conn.cursor() as cur:
		cur.execute('SELECT {f} FROM samples WHERE sample_id=\'{sample}\''.format(f=field, sample=params.sample))
		field = cur.fetchone()[0]
	return field

def zip_df(df):
	out = []
	for c in list(df.columns):
		out.append(df[c])	
	
	return zip(*out)

def get_compute_log(column, conn):
	with conn.cursor() as cur:
		cur.execute('SELECT {col} FROM compute_log WHERE sample_id=\'{sample}\''.format(col=column, sample=params.sample))
		field = cur.fetchone()[0]
	return field == 1

def compute_log_complete(job, params, conn):
	with conn.cursor() as cur:
		cur.execute('UPDATE compute_log SET {job} = %s WHERE sample_id=\'{sample}\''.format(job=job, sample=params.sample), (1,))

def insert_executemany(entries, table, conn):

	if len(entries) == 0:
		return None

	num_fields = len(entries[0])
	placeholders = '(' + ','.join(['%s']*num_fields) + ')'

	with conn.cursor() as cur:
		args_str = ','.join([cur.mogrify(placeholders, x).decode('UTF-8') for x in entries])
		cur.execute("INSERT INTO {tbl} VALUES ".format(tbl=table) + args_str)

def delete_where(table, column, equals, conn):

	with conn.cursor() as cur:
		cur.execute("DELETE FROM {tbl} WHERE {col}='{eq}'".format(tbl=table, col=column, eq=equals))

def fasta_gen(f):

	for rec in SeqIO.parse(f, 'fasta'):
		yield rec

def revcomp(seq):
	seq = seq.upper()
	rcomp = ''
	for c in seq[::-1]:
		if c == 'A':
			rcomp += 'T'
		elif c == 'T':
			rcomp += 'A'
		elif c == 'G':
			rcomp += 'C'
		elif c == 'C':
			rcomp += 'G'
		else:
			rcomp += c
	return rcomp

def canonical(seq):
	return sorted([seq.upper(), revcomp(seq)])[0]

def seqhash(seq, seqtype='nuc'):
	seq = seq.upper()

	if seqtype == 'nuc':	
		canon_seq = canonical(seq)
		return hashlib.sha256(canon_seq.encode('utf-8')).hexdigest()
	elif seqtype == 'prot':
		return hashlib.sha256(seq.strip('*').encode('utf-8')).hexdigest()

def gbucket_exists(f):
	storage_client = storage.Client()
	bucket = storage_client.bucket('durrant')
	blob = bucket.blob(f)
	return blob.exists()
	
def gbucket_download(f, o):
	storage_client = storage.Client()
	bucket = storage_client.bucket('durrant')
	blob = bucket.blob(f)
	blob.download_to_filename(o)

def gbucket_upload(f, o):
	storage_client = storage.Client()
	bucket = storage_client.bucket('durrant')
	blob = bucket.blob(o)
	blob.upload_from_filename(f)

def gunzip_file(fgz, outdir=None):
	f = fgz.replace('.gz', '')
	if outdir is not None:
		f = join(outdir, basename(f))
	outfile = open(f, 'w')
	with gzip.open(fgz, 'rt') as infile:
		for line in infile:
			outfile.write(line)
	outfile.close()
	return f

def gzip_file(f):
	fgz = f + '.gz'
	outfile = gzip.open(fgz, 'wt')
	with open(f) as infile:
		for line in infile:
			outfile.write(line)
	outfile.close()
	return fgz

def convert_gff_to_fasta(gff):
	if gff.endswith('.gz'):
		handle = gzip.open(gff, 'rt')
	else:
		handle = open(gff)

	with gzip.open(gff+'.tmp', 'wt') as outfile:
		fasta = False
		for line in handle:
			line = line.strip()
			if fasta:
				print(line, file=outfile)
			if line == '##FASTA':
				fasta = True
	handle.close()
	rename(gff+'.tmp', gff)

def download_genome(output, params):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		refseq_path = get_sample_field('refseq_path', params, conn)
		genbank_path = get_sample_field('genbank_path', params, conn)
		other_meta = get_sample_field('other_meta', params, conn)
		database = get_sample_field('database', params, conn)
		compute_log_done = get_compute_log('download_genome', conn)
		conn.close()

		download_dir = refseq_path
		if refseq_path == 'NaN':
			download_dir = genbank_path

		if database == 'uhgg':
			ftp_download_genome = other_meta.split('gff_download_path=')[-1].split(';')[-1]
		elif database in ['jgi', 'gem', 'tara', 'hgg', 'youngblut', 'mgrast']:
			ftp_download_genome = other_meta.split('dwnld_link=')[-1].split(';')[0]
			ftp_download_genome = ftp_download_genome.replace('gs://durrant/', '')
		elif database == 'mgnify':
			ftp_download_genome = other_meta.split('dwnld=')[-1].split(';')[0]
		else:
			ftp_download_genome = join(download_dir, basename(download_dir)+'_genomic.fna.gz')

		gbucket_download_genome = join("crispestdb/" + db_path, basename(db_path) + '.fna.gz')

		if gbucket_exists(gbucket_download_genome) and compute_log_done and params.rerun is False:
			print("Downloading existing genome...")
			gbucket_download(gbucket_download_genome, output.genome)
		else:
			print("Processing genome...")
			if database in ['jgi', 'gem', 'tara', 'hgg', 'youngblut', 'mgrast']:
				print("DOWNLOADING", ftp_download_genome)
				gbucket_download(ftp_download_genome, output.genome)
			else:
				command = 'wget --continue -t 3 -T 60 -w 5 -O {target} {source} &> /dev/null'.format(target=output.genome, source=ftp_download_genome)
				print('wget command:', command)
				shell(command)

			if database=='uhgg': convert_gff_to_fasta(output.genome)

			gbucket_upload(output.genome, gbucket_download_genome)

			conn = get_connection()
			compute_log_complete('download_genome', params, conn)
			conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if conn is not None:
			conn.close()

def pyfastx_genome(input, output):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('pyfastx_genome', conn)
		conn.close()

		gbucket_genome_pyfastx = join("crispestdb/" + db_path, basename(db_path), basename(output.genome_pyfastx))
		if gbucket_exists(gbucket_genome_pyfastx) and compute_log_done and params.rerun is False:
			print("Downloading existing genome index...")
			gbucket_download(gbucket_genome_pyfastx, output.genome_pyfastx)
			return

		cmd = 'pyfastx index {f}'.format(f=input.genome)
		print('command:', cmd)
		shell(cmd)

		gbucket_upload(output.genome_pyfastx, gbucket_genome_pyfastx)

		conn = get_connection()
		compute_log_complete('pyfastx_genome', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if conn is not None:
			conn.close()


def genome_stats(input, output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('genome_stats', conn)
		conn.close()
		
		gbucket_genome_stats = join("crispestdb/" + db_path, basename(output.genome_stats))
		gbucket_contigs = join("crispestdb/" + db_path, basename(output.contigs))

		if gbucket_exists(gbucket_genome_stats) and gbucket_exists(gbucket_contigs) and compute_log_done and params.rerun is False:
			print("Downloading existing genome_stats...")
			gbucket_download(gbucket_genome_stats, output.genome_stats)
			gbucket_download(gbucket_contigs, output.contigs)
		else:
			print("Processing genome_stats...")
			outfile = open(output.genome_stats, 'w')

			total_bases = 0
			contigs = []
			for rec in SeqIO.parse(gzip.open(input.genome, 'rt'), 'fasta'):
				contigs.append((rec.id, seqhash(str(rec.seq), seqtype='nuc'), rec.description, len(rec.seq)))
				total_bases += len(rec.seq)
			print('num_bases', total_bases, sep='\t', file=outfile)
			print('num_contigs', len(contigs), sep='\t', file=outfile)
			outfile.close()

			with gzip.open(output.contigs, 'wt') as outfile:
				print("contig_id", "seqhash", "description", "length", sep='\t', file=outfile)
				for line in contigs:
					print(*line, sep='\t', file=outfile)

			conn = get_connection()
			cur = conn.cursor()
			cur.execute('UPDATE samples SET num_bases = %s, num_contigs = %s WHERE sample_id=\'{sample}\''.format(sample=params.sample), (total_bases, len(contigs)))
			
			gbucket_upload(output.genome_stats, gbucket_genome_stats)
			gbucket_upload(output.contigs, gbucket_contigs)
			compute_log_complete('genome_stats', params, conn)
			conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if conn is not None:
			conn.close()


def rename_prodigal_output(faa, gff):
	faa_tmp = faa + '.tmp'
	out_faa_tmp = open(faa_tmp, 'w')

	name2seqhash = dict()
	for rec in SeqIO.parse(faa, 'fasta'):
		prot_seqhash = seqhash(str(rec.seq), seqtype='prot')
		name2seqhash[rec.id] = prot_seqhash
		rec.id = prot_seqhash
		rec.description = prot_seqhash
		SeqIO.write([rec], out_faa_tmp, 'fasta')

	out_faa_tmp.close()
	remove(faa)
	rename(faa_tmp, faa)
	
	
	gff_tmp = gff + '.tmp'
	out_gff_tmp = open(gff_tmp, 'w')
	with open(gff) as infile:
		for line in infile:
			if line.startswith('#'):
				out_gff_tmp.write(line)
				continue
			line = line.strip().split('\t')
			name = line[0] + '_' + line[-1].split(';')[0].split('_')[-1]

			line[-1] = re.sub(r'^ID=[^;]+', 'ID='+name2seqhash[name], line[-1])
			print('\t'.join(line), file=out_gff_tmp)

	out_gff_tmp.close()
	remove(gff)
	rename(gff_tmp, gff)

def count_records_in_fasta(fasta):
	records = 0
	if fasta.endswith('.gz'):
		with gzip.open(fasta, "rt") as infile:
			for line in infile:
				if line.startswith('>'):
					records += 1
	else:
		with open(fasta) as infile:
			for line in infile:
				if line.startswith('>'):
					records += 1
	return records

def combine_files(files, outfile, exclude_startswith=[]):
	with open(outfile, 'w') as out:

		for f in files:

			with open(f) as infile:
				for line in infile:

					skip = False
					for exclude in exclude_startswith:
						if line.startswith(exclude):
							skip = True
							break
					if skip is True:
						continue
					out.write(line)

	for f in files:
		remove(f)

def run_prodigal_simple(prodigal_path, infile, outprefix):
	FNULL = open(devnull, 'w')
	bashCommand = '{prodigal} -p meta -f gff -o {gff} -a {faa} -i {input}'.format(
		prodigal=prodigal_path, gff=outprefix + '.gff', faa=outprefix + '.faa',
		input=infile
	)

	print('prodigal command:', bashCommand)
	process = subprocess.Popen(bashCommand.split(), stdout=FNULL, stderr=FNULL)
	output, error = process.communicate()

def run_prodigal_multithread(prodigal_path, infile, outdir, gff_out, faa_out, threads):

	print("Counting records in FASTA file...")
	record_count = count_records_in_fasta(infile)
	print("The FASTA file contains %d records..." % record_count)

	print("Writing FASTA file to batches for multithreading..." )
	record_count_per_file = math.ceil(record_count / threads)
	filecount = 0
	outrecs = []
	outfiles = []

	if infile.endswith('.gz'):
		infile = gzip.open(infile, "rt")

	for record in SeqIO.parse(infile, 'fasta'):
		outrecs.append(record)
		if len(outrecs) == record_count_per_file:
			filecount += 1
			outfile = join(outdir, 'input%d.fna' % filecount)
			SeqIO.write(outrecs, outfile, 'fasta')
			outfiles.append(outfile)
			outrecs = []

	if len(outrecs) > 0:
		filecount += 1
		outfile = join(outdir, 'input%d.fna' % filecount)
		SeqIO.write(outrecs, outfile, 'fasta')
		outfiles.append(outfile)
	del outrecs

	prodigal_files = [join(outdir, 'prodigal%d' % (i+1)) for i in range(len(outfiles))]
	with Pool(processes=threads) as pool:
		pool.starmap(run_prodigal_simple, zip(cycle([prodigal_path]), outfiles, prodigal_files))

	combine_files(outfiles, join(outdir, 'input.fna'))
	combine_files([f+'.faa' for f in prodigal_files], faa_out)
	combine_files([f+'.gff' for f in prodigal_files], gff_out, ['#'])


def prodigal(input, output, params, threads):

	conn = None
	outdir = None
	try:

		conn = get_connection()
		domain = get_sample_field('domain', params, conn)
		num_bases = get_sample_field('num_bases', params, conn)
		db_path = get_sample_field('db_path', params, conn)
		sample_type = get_sample_field('sample_type', params, conn)
		compute_log_done = get_compute_log('prodigal', conn)
		conn.close()

		gbucket_faa = join("crispestdb/" + db_path, basename(db_path) + '.prodigal.faa.gz')
		gbucket_gff = join("crispestdb/" + db_path, basename(db_path) + '.prodigal.gff.gz')

		if not gbucket_exists(gbucket_faa) or not gbucket_exists(gbucket_gff) or not compute_log_done or params.rerun is True:

			outdir = output.faa + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)

			print("Running prodigal...")
			gunzipped_genome = gunzip_file(input.genome, outdir)
			tmp_gff = join(outdir, basename(output.gff.replace('.gz', '')))
			tmp_faa = join(outdir, basename(output.faa.replace('.gz', '')))

			if threads > 1 or sample_type in ['meta', 'metagenome'] or domain == 'viruses' or num_bases <= 20000:
				run_prodigal_multithread('prodigal', gunzipped_genome, outdir, tmp_gff, tmp_faa, threads)
			else:
				print("IMPLEMENT SINGLE THREAD PROKKA")
				command = 'prodigal -f gff -o {gff} -a {faa} -i {input} &> /dev/null'.format(
				    gff=tmp_gff, faa=tmp_faa, input=gunzipped_genome, 
				)
				print('prodigal command:', command)
				shell(command)

			rename_prodigal_output(tmp_faa, tmp_gff)
		
			final_faa, final_gff = gzip_file(tmp_faa), gzip_file(tmp_gff)

			rename(final_faa, output.faa)
			rename(final_gff, output.gff)
		
			gbucket_upload(output.faa, gbucket_faa)
			gbucket_upload(output.gff, gbucket_gff)

			gene_count = 0
			for rec in SeqIO.parse(gzip.open(output.faa, 'rt'), "fasta"):
				gene_count += 1

			conn = get_connection()
			with conn.cursor() as cur:
				cur.execute('UPDATE samples SET num_genes_prodigal = %s WHERE sample_id=\'{sample}\''.format(sample=params.sample), (gene_count,))

			compute_log_complete('prodigal', params, conn)
			conn.commit()

		else:
			print("Downloading existing prodigal files...")
			gbucket_download(gbucket_faa, output.faa)
			gbucket_download(gbucket_gff, output.gff)


	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def hmmsearch_crisprcastyper(input, output, params, threads):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('hmmsearch_crisprcastyper', conn)
		conn.close()	

		gbucket_hmmsearch_crisprcastyper = join("crispestdb/" + db_path, basename(output.hmmsearch_out))

		if gbucket_exists(gbucket_hmmsearch_crisprcastyper) and compute_log_done and params.rerun is False:
			print("Downloading existing hmmsearch CRISPRCasTyper files...")
			gbucket_download(gbucket_hmmsearch_crisprcastyper, output.hmmsearch_out)
		else:
			outdir = output.hmmsearch_out + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)
			makedirs('data/crisprcastyper_profiles.hmm')
			if not isfile('data/crisprcastyper_profiles.hmm'):
				gbucket_crisprcastyper_profiles_hmm = join("crispestdb-other", 'crisprcastyper_profiles.hmm')
				gbucket_download(gbucket_crisprcastyper_profiles_hmm, 'data/crisprcastyper_profiles.hmm')
				
			gunzipped_faa = gunzip_file(input.faa, outdir)
			outtbl = join(outdir, basename(output.hmmsearch_out.replace('.tsv.gz', '')+'.tbl'))
			print("Running hmmsearch...")
			command = 'hmmsearch --cpu {threads} -E 1 --tblout {outtbl} {hmmfile} {seqdb} &> /dev/null'.format(
				threads=threads, outtbl=outtbl, hmmfile='data/crisprcastyper_profiles.hmm', seqdb=gunzipped_faa
			)
			print('hmmsearch command:', command)
			shell(command)
	
			outfile = gzip.open(output.hmmsearch_out, 'wt')
			print('target_name', 'query_name', 'full_seq_evalue', 'full_seq_score', 'full_seq_bias', 'best_domain_evalue', 'best_domain_score', 'best_domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', sep='\t', file=outfile)
			keep_genes = set()
			hmmsearch_to_database = []
			header = ['target_name', 'target_accession', 'query_name', 'query_accession', 'full_seq_evalue', 'full_seq_score', 'full_seq_bias', 'best_domain_evalue', 'best_domain_score', 'best_domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'target_description']
			with open(outtbl) as infile:
				for line in infile:
					if line.startswith('#'):
						continue
					line = line.strip().split()
					line = {header[i]:line[i] for i in range(len(header))}
					keep_genes.add(line['target_name'])
					hmmsearch_to_database.append((params.sample, line['target_name'], line['query_name'], float(line['full_seq_evalue']), float(line['full_seq_score']), float(line['best_domain_evalue']), float(line['best_domain_score'])))
					print(line['target_name'], line['query_name'], line['full_seq_evalue'], line['full_seq_score'], line['full_seq_bias'], line['best_domain_evalue'], line['best_domain_score'], line['best_domain_bias'], line['exp'], line['reg'], line['clu'], line['ov'], line['env'], line['dom'], line['rep'], line['inc'], sep='\t', file=outfile)
			outfile.close()
	
			gbucket_upload(output.hmmsearch_out, gbucket_hmmsearch_crisprcastyper)

			conn = get_connection()
			compute_log_complete('hmmsearch_crisprcastyper', params, conn)
			conn.commit()
	
	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def mash_sketch(input, output, params, threads):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('mash_sketch', conn)
		conn.close()

		gbucket_mash = join("crispestdb/" + db_path, basename(output.mash_out))

		if gbucket_exists(gbucket_mash) and compute_log_done and params.rerun is False:
			print("Downloading existing mash file...")
			gbucket_download(gbucket_mash, output.mash_out)
		else:
			print("Running mash index...")

			outdir = output.mash_out + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)
	
			gunzipped_genome = gunzip_file(input.genome, outdir)
			command = 'mash sketch -p {threads} {genome} &> /dev/null'.format(threads=threads, genome=gunzipped_genome)
			print('mash command:', command)
			shell(command)
			rename(gunzipped_genome+'.msh', output.mash_out)
	
			gbucket_upload(output.mash_out, gbucket_mash)

			conn = get_connection()
			compute_log_complete('mash_sketch', params, conn)
			conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def genomesearch_markers(input, output, params, threads):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('genomesearch_markers', conn)
		conn.close()
	
		gbucket_genomesearch_markers = join("crispestdb/" + db_path, basename(output.genomesearch_markers))
	
		if gbucket_exists(gbucket_genomesearch_markers) and compute_log_done and params.rerun is False:
			print("Downloading existing genomesearch markers...")
			gbucket_download(gbucket_genomesearch_markers, output.genomesearch_markers)
		else:
			print("Finding genomesearch marker genes...")
	
			outdir = output.genomesearch_markers + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)
			makedirs('data', exist_ok=True)
			if not isfile('data/phylophlan_marker_references.dmnd'):
				gbucket_phylophlan_marker_references = join("crispestdb-other", 'phylophlan_marker_references.dmnd')
				gbucket_download(gbucket_phylophlan_marker_references, 'data/phylophlan_marker_references.dmnd')

			diamond_results = join(outdir, basename(output.genomesearch_markers+'.dmd.tsv'))
			command = 'diamond blastp --query {0} --out {1} --outfmt 6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore --db {2} --threads {3} &> /dev/null'.format(input.faa, diamond_results, 'data/phylophlan_marker_references.dmnd', threads)
			print('diamond command:', command)
			shell(command)
	
			top_markers = dict()
			with open(diamond_results) as infile:
				for line in infile:
					qseqid, sseqid, qlen, slen, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore =  line.strip().split('\t')
					qlen, slen, qstart, qend, sstart, send = int(qlen), int(slen), int(qstart), int(qend), int(sstart), int(send)
					length, evalue, bitscore = int(length), float(evalue), float(bitscore)
					finding = (qseqid, length, evalue, bitscore)
					marker = sseqid.split('_')[1]
	
					qaln = qend - qstart
					saln = send - sstart
	
					query_smaller = True
					if slen > qlen:
						query_smaller = False
	
					if finding[-2] >= 1e-6:
						continue
	
					if min([qlen, slen]) / max([qlen, slen]) < 0.85:
						continue
	
					if query_smaller:
						alnlength = qaln
					else:
						alnlength = saln
	
					if alnlength / min([qlen, slen]) < 0.85:
						continue
	
					if marker not in top_markers:
						top_markers[marker] = finding
					else:
						if finding[-1] > top_markers[marker][-1]:
							top_markers[marker] = finding
	
			marker2gene = dict()
			gene2marker = dict()
			for rec in top_markers:
				marker2gene[rec] = top_markers[rec][0]
				gene2marker[top_markers[rec][0]] = rec
	
			outfile = gzip.open(output.genomesearch_markers, 'wt')
			print("gene", "marker", sep='\t', file=outfile)
			for gene in gene2marker:
				print(gene, gene2marker[gene], sep='\t', file=outfile)
			outfile.close()
	
			gbucket_upload(output.genomesearch_markers, gbucket_genomesearch_markers)

			conn = get_connection()
			compute_log_complete('genomesearch_markers', params, conn)
			conn.commit()
	
	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def run_minced_simple(infasta, outdir):

	inprefix = join(dirname(infasta), basename(infasta).split('.')[0])
	tmp_minced_out = inprefix + '.out'
	tmp_minced_gff = inprefix + '.gff'

	command = 'minced -spacers -minNR 1 -minSL 20 -gffFull {genome} {out} {gff}'.format(genome=infasta, out=tmp_minced_out, gff=tmp_minced_gff)
	print('minced command:', command)
	shell(command)

	out_gff = []
	with open(tmp_minced_gff) as infile:
		for line in infile:
			if line.startswith("#"):
				continue
			out_gff.append(line.strip().split('\t'))

	tmp_repeat_units_gff = join(outdir, basename(tmp_minced_gff + '.repeat_units.gff'))
	tmp_repeat_units_fna = join(outdir, basename(tmp_minced_gff + '.repeat_units.fna'))
	with open(tmp_repeat_units_gff, 'w') as outfile:
		for line in out_gff:
			if line[2] == 'repeat_unit':
				print(*line, sep='\t', file=outfile)

	command = 'bedtools getfasta -fi {fasta} -bed {units} -fo {units_fna} &> /dev/null'.format(units=tmp_repeat_units_gff, fasta=infasta, units_fna=tmp_repeat_units_fna)
	print('command:', command)
	shell(command)

	repeat_units = fasta_gen(tmp_repeat_units_fna)
	spacers = fasta_gen(inprefix + '_spacers.fa')

	tmp_minced_out_final = inprefix + '.minced.gff'
	outfile = open(tmp_minced_out_final, 'w')
	for i, line in enumerate(out_gff):

		group = line[-1].split(';')[0].split('=')[-1]
		if line[2] == 'repeat_region':
			seq = line[-1].split('=')[-1]
		elif line[2] == 'repeat_unit':
			seq = str(next(repeat_units).seq)
			line[-1] = line[-1] + ';seq=' + seq
		
		print(*line, sep='\t', file=outfile)

		if line[2] == 'repeat_unit' and i+1 != len(out_gff) and out_gff[i+1][2] != 'repeat_region':

			spacer_rec = next(spacers)
			seq = str(spacer_rec.seq)

			contig, start, end = line[0], str(int(line[4])+1), str(int(out_gff[i+1][3])-1)

			line[0], line[2], line[3], line[4] = contig, 'spacer', start, end
			line[-1] = line[-1].split(';seq=')[0] + ';seq=' + seq
			print(*line, sep='\t', file=outfile)

	outfile.close()

	return tmp_minced_out_final


def run_minced_multithread(infile, final_outfile, outdir, threads):

	print("Counting records in FASTA file...")
	record_count = count_records_in_fasta(infile)
	print("The FASTA file contains %d records..." % record_count)

	print("Writing FASTA file to batches for multithreading..." )
	record_count_per_file = math.ceil(record_count / threads)
	filecount = 0
	outrecs = []
	outfiles = []

	if infile.endswith('.gz'):
		infile = gzip.open(infile, "rt")

	for record in SeqIO.parse(infile, 'fasta'):
		outrecs.append(record)
		if len(outrecs) == record_count_per_file:
			filecount += 1
			outfile = join(outdir, 'input%d.fna' % filecount)
			SeqIO.write(outrecs, outfile, 'fasta')
			outfiles.append(outfile)
			outrecs = []

	if len(outrecs) > 0:
		filecount += 1
		outfile = join(outdir, 'input%d.fna' % filecount)
		SeqIO.write(outrecs, outfile, 'fasta')
		outfiles.append(outfile)
	del outrecs

	args = [(outf, outdir) for outf in outfiles]
	with Pool(processes=16) as pool:
		outfiles = pool.starmap(run_minced_simple, args)

	crispr_count = 0
	with gzip.open(final_outfile, 'wt') as outfile:
		for f in outfiles:
			crispr_count += 1
			corrected_crispr_name = 'CRISPR' + str(crispr_count)
			current_crispr_name = 'CRISPR1'
			counter = 0
			with open(f) as infile:
				for line in infile:
					counter += 1
					next_crispr_name = re.findall('CRISPR\d+|$', line)[0]
					if current_crispr_name != next_crispr_name:
						crispr_count += 1
						corrected_crispr_name = 'CRISPR' + str(crispr_count)
					current_crispr_name = next_crispr_name
					line = re.sub('CRISPR\d+', corrected_crispr_name, line)
					outfile.write(line)

			if counter == 0:
				crispr_count -= 1
				corrected_crispr_name = 'CRISPR' + str(crispr_count)


def minced(input, output, params, threads):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('minced', conn)
		conn.close()
	
		gbucket_minced = join("crispestdb/" + db_path, basename(output.minced))
		if gbucket_exists(gbucket_minced) and compute_log_done and params.rerun is False:
			print("Downloading existing minced file...")
			gbucket_download(gbucket_minced, output.minced)
		else:
			print("Running minced...")
			
			outdir = output.minced + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)
	
			run_minced_multithread(input.genome, output.minced, outdir, threads)
	
			gbucket_upload(output.minced, gbucket_minced)

			conn = get_connection()
			compute_log_complete('minced', params, conn)
			conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def run_crisprcasfinder_simple(infasta, outdir, add_pwd=True):

	inprefix = join(dirname(infasta), basename(infasta).split('.')[0])
	tmp_outdir = inprefix + '.outdir'
	makedirs(tmp_outdir, exist_ok=True)

	if add_pwd:
		command = 'cd {outdir}; /usr/bin/perl /home/mdurrant/CRISPRCasFinder/CRISPRCasFinder.pl -so /home/mdurrant/CRISPRCasFinder/sel392v2.so -cf /home/mdurrant/CRISPRCasFinder/CasFinder-2.0.3 -drpt /home/mdurrant/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /home/mdurrant/CRISPRCasFinder/supplementary_files/Repeat_List.csv -def G -in {pwd}/{genome} -minSP 20'.format(pwd=pwd, outdir=tmp_outdir, genome=infasta)
	else:
		command = 'cd {outdir}; /usr/bin/perl /home/mdurrant/CRISPRCasFinder/CRISPRCasFinder.pl -so /home/mdurrant/CRISPRCasFinder/sel392v2.so -cf /home/mdurrant/CRISPRCasFinder/CasFinder-2.0.3 -drpt /home/mdurrant/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /home/mdurrant/CRISPRCasFinder/supplementary_files/Repeat_List.csv -def G -in {genome} -minSP 20'.format(outdir=tmp_outdir, genome=infasta)
		
	print('crisprcasfinder command:', command)
	shell(command)

	out_gff = []
	for f in glob(join(tmp_outdir, '*', 'GFF', '*gff')):
		with open(f) as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				line = line.strip().split('\t')	
				if len(line) == 1:
					continue
				if line[2] == 'CRISPRspacer':
					line[2] = 'spacer'
				elif line[2] == 'CRISPRdr':
					line[2] = 'repeat_unit'
				group = line[-1].split('Parent=')[-1].split(';')[0]
				out_gff.append(line)
			
				seq = line[-1].split(';')[0].split('=')[-1]

	outfile_path = inprefix + '.crisprcasfinder.gff'
	outfile = open(outfile_path, 'w')
	for i, line in enumerate(out_gff):
		print(*line, sep='\t', file=outfile)
	outfile.close()

	return outfile_path


def run_crisprcasfinder_multithread(infile, final_outfile, outdir, threads):

	print("Counting records in FASTA file...")
	record_count = count_records_in_fasta(infile)
	print("The FASTA file contains %d records..." % record_count)

	print("Writing FASTA file to batches for multithreading..." )
	record_count_per_file = math.ceil(record_count / threads)
	filecount = 0
	outrecs = []
	outfiles = []

	if infile.endswith('.gz'):
		infile = gzip.open(infile, "rt")

	for record in SeqIO.parse(infile, 'fasta'):
		outrecs.append(record)
		if len(outrecs) == record_count_per_file:
			filecount += 1
			outfile = join(outdir, 'input%d.fna' % filecount)
			SeqIO.write(outrecs, outfile, 'fasta')
			outfiles.append(outfile)
			outrecs = []

	if len(outrecs) > 0:
		filecount += 1
		outfile = join(outdir, 'input%d.fna' % filecount)
		SeqIO.write(outrecs, outfile, 'fasta')
		outfiles.append(outfile)
	del outrecs

	args = [(outf, outdir) for outf in outfiles]
	with Pool(processes=threads) as pool:
		outfiles = pool.starmap(run_crisprcasfinder_simple, args)

	with gzip.open(final_outfile, 'wt') as outfile:
		for f in outfiles:
			with open(f) as infile:
				for line in infile:
					outfile.write(line)

def crisprcasfinder(input, output, params, threads):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('crisprcasfinder', conn)
		conn.close()
	
		gbucket_crisprcasfinder = join("crispestdb/" + db_path, basename(output.crisprcasfinder))
		if gbucket_exists(gbucket_crisprcasfinder) and compute_log_done and params.rerun is False:
			print("Downloading existing CRISPRCasFinder file...")
			gbucket_download(gbucket_crisprcasfinder, output.crisprcasfinder)
	
		else:
			print("Running CRISPRCasFinder...")
			outdir = output.crisprcasfinder + '.outdir'
			if os.path.isdir(outdir):
				shutil.rmtree(outdir)
			makedirs(outdir, exist_ok=True)
			
			run_crisprcasfinder_multithread(input.genome, output.crisprcasfinder, outdir, threads)
	
			gbucket_upload(output.crisprcasfinder, gbucket_crisprcasfinder)

			conn = get_connection()
			compute_log_complete('crisprcasfinder', params, conn)
			conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def find_cas_fusions(input, output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('find_cas_fusions', conn)
		conn.close()
	
		gbucket_casfusions_summary = join("crispestdb/" + db_path, basename(output.cas_fusions_bed))
		gbucket_casfusions_bed = join("crispestdb/" + db_path, basename(output.cas_fusions_summary))
		gbucket_crisprcastyper_domains = join("crispestdb/" + db_path, basename(output.crisprcastyper_domains))
		if gbucket_exists(gbucket_casfusions_summary) and gbucket_exists(gbucket_casfusions_bed) and gbucket_exists(gbucket_crisprcastyper_domains) and compute_log_done and params.rerun is False:
			print("Downloading existing Cas Fusions files...")
			gbucket_download(gbucket_casfusions_bed, output.cas_fusions_bed)
			gbucket_download(gbucket_casfusions_summary, output.cas_fusions_summary)
			gbucket_download(gbucket_crisprcastyper_domains, output.crisprcastyper_domains)
			return	

		outdir = output.cas_fusions_bed + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)

		crispr_genes = defaultdict(set)
		#gene_int = ['RT_0_CAS-I-II-III-IV-V-VI', 'RT_0_CAS-I-E', 'RT_0_CAS-III-AM-III-B', 'gRAMP_0_CAS-III-E']

		#filter and save crispr hits
		with gzip.open(input.hmmsearch_out, 'rt') as file1:
			reader1 = csv.reader(file1, delimiter='\t')
			#skip header
			next(reader1)
			for line in reader1:
				#filter e-value
				if float(line[2]) > 1e-4:
					continue
				full_query = line[1]
				query = full_query.split('|')[1]
				#if 'Cas' in query or query in gene_int:
				head = line[0].rstrip()
				crispr_genes[head].add(full_query)

		#save fasta sequences somewhere
		with gzip.open(input.faa, 'rt') as file2:
			for rec in SeqIO.parse(file2, 'fasta'):
				if rec.id not in crispr_genes:
					continue
				for name in crispr_genes[rec.id]:
					name = name.replace('|', '_')
					outfile = open(join(outdir, name+'.filtered.faa'), 'a')
					SeqIO.write([rec], outfile, 'fasta')
					outfile.close()

		#hmmsearch
		start = timeit.default_timer()
		for f in glob(join(outdir, '*.filtered.faa')):
			hmm = basename(f).split('.')[0]
			hmm_file = join('/home/mdurrant/data/crisprcastyper_profiles', hmm+'.hmm')
			hmm_out = join(outdir, hmm+'.hmmsearch.tbl')
			shell('hmmsearch --domtblout %s %s %s &> /dev/null' % (hmm_out, hmm_file, f))
		stop = timeit.default_timer()
		print('total run time:', (stop - start), 'seconds')

		domain_header = ["target_name", "tacc", "tlen", "query_name", "qacc", "qlen", "full_seq_evalue", "full_seq_score", "full_seq_bias", "domain_number", "domain_of", "domain_c_evalue", "domain_i_evalue", "domain_score", "domain_bias", "hmm_from", "hmm_to", "align_from", "align_to", "env_from", "env_to", "acc", "desc"]
		outfile = gzip.open(output.crisprcastyper_domains, 'wt')
		print(*domain_header, sep='\t', file=outfile)
		for f in glob(join(outdir, '*.hmmsearch.tbl')):
			with open(f) as infile:
				for line in infile:
					if line.startswith("#"):
						continue
					line = line.strip().split()
					print(*line, sep='\t', file=outfile)
		outfile.close()

		bedfile = join(outdir, 'crisprcastyper_domains.bed')
		outfile = open(bedfile, 'w')
		with gzip.open(output.crisprcastyper_domains, 'rt') as infile:
			header = infile.readline().strip().split()
			for line in infile:
				line = line.strip().split()
				line = {header[i]:line[i] for i in range(len(line))}
				Name = 'Name='+line['query_name']+';cevalue='+line['domain_c_evalue']+';dscore='+line['domain_score']+';hmm_from='+line['hmm_from']+';hmm_to='+line['hmm_to']
				print(line['target_name'], line['align_from'], line['align_to'], Name, sep='\t', file=outfile)
		outfile.close()

		intersect_bedfile = join(outdir, 'crisprcastyper_domains.intersect.bed')
		shell("bedtools intersect -a {bed} -b {bed} -wa -wb -f 0.9 -r | awk '$4 != $8' > {outbed}".format(bed=bedfile, outbed=intersect_bedfile))

		try:
			df = pd.read_csv(intersect_bedfile, sep='\t', header=None)
			df.columns = ['ref1', 'start1', 'end1', 'name1', 'ref2', 'start2', 'end2', 'name2']
		except EmptyDataError:
			df = pd.DataFrame(columns=['ref1', 'start1', 'end1', 'name1', 'ref2', 'start2', 'end2', 'name2'])

		unique_entries = set([(ref1, start1, end1, name1) for ref1, start1, end1, name1 in zip(df.ref1, df.start1, df.end1, df.name1)])
		exclude_entries = set()
		keep_entries = set()
		for ref1, start1, end1, name1 in unique_entries:
			score1 = float(name1.split('dscore=')[-1].split(';')[0])
			df_filt = df[(df.ref1 == ref1) & (df.start1 == start1) & (df.end1 == end1) & (df.name1 == name1)]
			keep = True
			equals = set()
			for ref2, start2, end2, name2 in zip(df_filt.ref2, df_filt.start2, df_filt.end2, df_filt.name2):
				score2 = float(name2.split('dscore=')[-1].split(';')[0])
				if score2 > score1:
					keep = False
				if score2 == score1:
					equals.add((ref2, start2, end2, name2))
				
			if keep:
				for entry in equals:
					if entry not in keep_entries:
						exclude_entries.add(entry)
				keep_entries.add((ref1, start1, end1, name1))
				
		try:
			df = pd.read_csv(bedfile, sep='\t', header=None)
			df.columns = ['ref', 'start', 'end', 'name']
		except:
			df = pd.DataFrame(columns=['ref', 'start', 'end', 'name'])
			
		filtered_bedfile = join(outdir, 'crisprcastyper_domains.filtered.bed')
		outbed = open(filtered_bedfile, 'w')
		for entry in zip(df.ref, df.start, df.end, df.name):
			if entry not in unique_entries:
				print(*entry, sep='\t', file=outbed)
				continue
			if entry in keep_entries and entry not in exclude_entries:
				print(*entry, sep='\t', file=outbed)
		outbed.close()

		shell('bedtools sort -i {bed} | bedtools merge | bedtools intersect -a stdin -b {bed} -wa -wb | gzip > {outbed}'.format(bed=filtered_bedfile, outbed=output.cas_fusions_bed))
		
		region2hmm = defaultdict(lambda: defaultdict(list))
		with gzip.open(output.cas_fusions_bed, 'rt') as infile:
	
			for line in infile:
				line = line.strip().split()
				name = line[-1]
				region = line[0]+':'+line[1]+'-'+line[2]
				region2hmm[line[0]][region].append(name.split('Name=')[-1].split(';')[0] + '(' + line[4] + '-' + line[5] + ')')

		outfile = gzip.open(output.cas_fusions_summary, 'wt')
		for name in region2hmm:
			out = ''
			for region in region2hmm[name]:
				out += region+'('+';'.join(region2hmm[name][region])+');'

			print(name, out, sep='\t', file=outfile)
		outfile.close()

		gbucket_upload(output.cas_fusions_bed, gbucket_casfusions_bed)
		gbucket_upload(output.cas_fusions_summary, gbucket_casfusions_summary)
		gbucket_upload(output.crisprcastyper_domains, gbucket_crisprcastyper_domains)

		conn = get_connection()
		compute_log_complete('find_cas_fusions', params, conn)
		conn.commit()
					
	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def hmmsearch_isescan(input, output, params, threads):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('hmmsearch_isescan', conn)
		conn.close()	

		gbucket_hmmsearch_isescan = join("crispestdb/" + db_path, basename(output.hmmsearch_isescan_out))

		if gbucket_exists(gbucket_hmmsearch_isescan) and compute_log_done and params.rerun is False:
			print("Downloading existing hmmsearch ISEScan files...")
			gbucket_download(gbucket_hmmsearch_isescan, output.hmmsearch_isescan_out)
			return

		outdir = output.hmmsearch_isescan_out + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)

		gunzipped_faa = gunzip_file(input.faa, outdir)
		outtbl = join(outdir, basename(output.hmmsearch_isescan_out.replace('.tsv.gz', '')+'.tbl'))
		print("Running hmmsearch...")
		command = 'hmmsearch --cpu {threads} -E 1 --domtblout {outtbl} {hmmfile} {seqdb} &> /dev/null'.format(
			threads=threads, outtbl=outtbl, hmmfile=params.isescan_hmm_path, seqdb=gunzipped_faa
		)
		print('hmmsearch command:', command)
		shell(command)

		outfile = gzip.open(output.hmmsearch_isescan_out, 'wt')
		domain_header = ["target_name", "tacc", "tlen", "query_name", "qacc", "qlen", "full_seq_evalue", "full_seq_score", "full_seq_bias", "domain_number", "domain_of", "domain_c_evalue", "domain_i_evalue", "domain_score", "domain_bias", "hmm_from", "hmm_to", "align_from", "align_to", "env_from", "env_to", "acc", "desc"]
		print(*domain_header, sep='\t', file=outfile)
		keep_genes = set()
		hmmsearch_to_database = []
		header = ['target_name', 'target_accession', 'query_name', 'query_accession', 'full_seq_evalue', 'full_seq_score', 'full_seq_bias', 'best_domain_evalue', 'best_domain_score', 'best_domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'target_description']
		with open(outtbl) as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				line = line.strip().split()
				print(*line, sep='\t', file=outfile)
		outfile.close()

		gbucket_upload(output.hmmsearch_isescan_out, gbucket_hmmsearch_isescan)

		conn = get_connection()
		compute_log_complete('hmmsearch_isescan', params, conn)
		conn.commit()
	
	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def makeblastdb(input, output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('makeblastdb', conn)
		conn.close()	

		gbucket_nhr = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nhr))
		gbucket_nin = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nin))
		gbucket_nog = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nog))
		gbucket_nsd = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nsd))
		gbucket_nsi = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nsi))
		gbucket_nsq = join("crispestdb/" + db_path, 'blastdb', basename(output.blastdb_nsq))

		if gbucket_exists(gbucket_nhr) and gbucket_exists(gbucket_nsq) and compute_log_done and params.rerun is False:
			print("Downloading existing blast db...")
			gbucket_download(gbucket_nhr, output.blastdb_nhr)
			gbucket_download(gbucket_nin, output.blastdb_nin)
			gbucket_download(gbucket_nog, output.blastdb_nog)
			gbucket_download(gbucket_nsd, output.blastdb_nsd)
			gbucket_download(gbucket_nsi, output.blastdb_nsi)
			gbucket_download(gbucket_nsq, output.blastdb_nsq)
			return

		shell('gunzip -k -f {genome}'.format(genome=input.genome))
		genome_nogz = input.genome.replace('.gz', '')

		blastdb_pref = output.blastdb_nhr.replace('.nhr', '')
		shell('makeblastdb -dbtype nucl -in {genome} -parse_seqids -out {blastdb}'.format(genome=genome_nogz, blastdb=blastdb_pref))

		gbucket_upload(output.blastdb_nhr, gbucket_nhr)
		gbucket_upload(output.blastdb_nin, gbucket_nin)
		gbucket_upload(output.blastdb_nog, gbucket_nog)
		gbucket_upload(output.blastdb_nsd, gbucket_nsd)
		gbucket_upload(output.blastdb_nsi, gbucket_nsi)
		gbucket_upload(output.blastdb_nsq, gbucket_nsq)

		remove(genome_nogz)

		conn = get_connection()
		compute_log_complete('makeblastdb', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if conn is not None:
			conn.close()

def self_targeting_spacers(input, output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('self_targeting_spacers', conn)
		conn.close()	

		gbucket_sts = join("crispestdb/" + db_path, basename(output.self_targeting_spacers))

		if gbucket_exists(gbucket_sts) and compute_log_done and params.rerun is False:
			print("Downloading existing self-targeting spacers...")
			gbucket_download(gbucket_sts, output.self_targeting_spacers)
			return

		outdir = output.self_targeting_spacers + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)

		shell('gunzip -k -f {genome}'.format(genome=input.genome))
		genome_nogz = input.genome.replace('.gz', '')
		
		tmp_fasta = join(outdir, 'tmp.fna')
		shell("zcat {gff} | awk '$3 == \"spacer\"' | bedtools getfasta -fi {genome} -bed stdin > {outfile}".format(genome=genome_nogz, gff=input.minced, outfile=tmp_fasta))
		shell("zcat {gff} | awk '$3 == \"spacer\"' | bedtools getfasta -fi {genome} -bed stdin >> {outfile}".format(genome=genome_nogz, gff=input.crisprcasfinder, outfile=tmp_fasta))

		outrecs = dict()
		for rec in SeqIO.parse(tmp_fasta, 'fasta'):
			outrecs[rec.id] = rec

		spacers_fasta = join(outdir, 'spacers.fna')
		SeqIO.write(list(outrecs.values()), spacers_fasta, 'fasta')
		
		workdir = dirname(output.self_targeting_spacers)
		blastprefix = join(workdir, params.sample+'.blastdb')
		outblast = join(outdir, 'blast.tsv')
		shell("blastn -task blastn-short -evalue 1 -query {spacers_fna} -db {blast_pref} -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident evalue bitscore qseq sseq' > {out}".format(spacers_fna=spacers_fasta, out=outblast, blast_pref=blastprefix))


		df_header = ['qseqid', 'sseqid', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'pident', 'evalue', 'bitscore', 'qseq', 'sseq']
		try:
			df = pd.read_csv(outblast, sep='\t', header=None)
			df.columns = df_header
		except EmptyDataError:
			df = pd.DataFrame(columns=df_header)

		try:
			df['sseqid'] = [s.split('|')[1] for s in df['sseqid']]
		except:
			pass
	

		outsts = join(outdir, 'sts.tsv')
		outfile = open(outsts, 'w')
		for qseqid, sseqid, length, qlen, slen, qstart, qend, sstart, send, pident, evalue, bitscore, qseq, sseq in zip_df(df):
			qseqid_orig, qstart_orig, qend_orig = qseqid.split(':')[0], int(qseqid.split(':')[-1].split('-')[0]), int(qseqid.split(':')[-1].split('-')[1])
			
			if qstart_orig+1 == sstart and qend_orig == send:
				continue

			if sstart < send:
				print(sseqid, sstart-1, send, "Name="+qseqid+";qlen="+str(qlen)+";qstart="+str(qstart)+";qend="+str(qend)+";pident="+str(pident)+";evalue="+str(evalue)+';qseq='+qseq+';sseq='+sseq, sep='\t', file=outfile)
			else:
				print(sseqid, send-1, sstart, "Name="+qseqid+";qlen="+str(qlen)+";qstart="+str(qstart)+";qend="+str(qend)+";pident="+str(pident)+";evalue="+str(evalue)+';qseq='+qseq+';sseq='+sseq, sep='\t', file=outfile)
		outfile.close()

		gff_sorted = join(outdir, 'sorted.gff')
		shell('bedtools sort -i {gff} > {outfile}'.format(gff=input.gff, outfile=gff_sorted))
		out_closest = join(outdir, 'closest.bed')
		shell('bedtools sort -i {infile} | bedtools closest -k 2 -D b -t all -a stdin -b {gff} > {out}'.format(infile=outsts, gff=gff_sorted, out=out_closest))

		outfile = gzip.open(output.self_targeting_spacers, 'wt')
		print("contig", "start", "end", "qseqid", "qlen", "qstart", "qend", "pident", "evalue", "qseq_aln", "sseq_aln", "dist", "gene_loc", "gene", sep='\t', file=outfile)
		with open(out_closest) as infile:
			for line in infile:
				line = line.strip().split('\t')
				gene_loc = line[4] + ':'+ line[7]+'-'+line[8]
				protein_id = line[-2].replace('ID=', '').split(';')[0]
				name = dict([l.split('=') for l in line[3].split(';')])
				print(line[0], line[1], line[2], name['Name'], name['qlen'], name['qstart'], name['qend'], name['pident'], name['evalue'], name['qseq'], name['sseq'], line[-1], gene_loc, protein_id, sep='\t', file=outfile)
		outfile.close()
		
		gbucket_upload(output.self_targeting_spacers, gbucket_sts)
		
		conn = get_connection()
		compute_log_complete('self_targeting_spacers', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def query_contig_length(contig):
	cmd = 'SELECT length FROM contig_lengths WHERE contig="{c}"'.format(c=contig)
	clength = cur.execute(cmd).fetchone()[0]
	return clength
	
def memoize(f):
	global contig_length_memo
	contig_length_memo = {}
	def helper(x):
		if x not in contig_length_memo:            
			contig_length_memo[x] = f(x)
		return contig_length_memo[x]
	return helper

def make_cas_neighbors_df(cas_intersect, cas_genes, sqldb):

	conn = sqlite3.connect(sqldb)
	global cur
	cur = conn.cursor()
	contig_length = memoize(query_contig_length)
	out = list()
	with open(cas_intersect) as infile:
		for line in infile:
			line = line.strip().split('\t')
			contig, bait_start, bait_end, bait_orient = line[0], int(line[3])-1, int(line[4]), line[5]
			bait_id = line[6].split('ID=')[-1].split(';')[0]
			target_start, target_end, target_orient = int(line[10])-1, int(line[11]), line[13]
			target_id = line[-1].split('ID=')[-1].split(';')[0]

			if bait_id == target_id and bait_start == target_start and bait_end == target_end:
				continue

			dist = -sys.maxsize
			if target_start < bait_end and bait_start < target_end:
				dist = 0
			elif target_orient == '+' and bait_start >= target_end:
				dist = bait_start - target_end + 1
			elif target_orient == '+' and bait_end <= target_start:
				dist = (target_start - bait_end + 1)*-1
			elif target_orient == '-' and bait_start >= target_end:
				dist = (bait_start - target_end + 1)*-1
			elif target_orient == '-' and bait_end <= target_start:
				dist = target_start - bait_end + 1
			elif bait_start >= target_end:
				dist = bait_start - target_end + 1
			elif bait_end <= target_start:
				dist = target_start - bait_end + 1
				
			clen = contig_length(contig)

			out.append([target_id[:24], bait_id[:24], contig, clen, target_start, target_end, target_orient, bait_start, bait_end, bait_orient, dist, 'Cas', cas_genes[bait_id][0], cas_genes[bait_id][1], cas_genes[bait_id][2]])

	cur.close()
	conn.close()
	df = pd.DataFrame(out)
	colnames = ['target_id', 'bait_id', 'contig', 'contig_length', 'target_start', 'target_end', 'target_orient', 'bait_start', 'bait_end', 'bait_orient', 'dist', 'bait_type', 'bait_domain', 'bait_evalue', 'bait_score']
	try:
		df.columns = colnames
	except:
		df = pd.DataFrame(columns=colnames)
	return df

def make_crispr_neighbors_df(crispr_intersect, crispr_info, sqldb):

	conn = sqlite3.connect(sqldb)
	global cur
	cur = conn.cursor()
	contig_length = memoize(query_contig_length)
	out = list()
	with open(crispr_intersect) as infile:
		for line in infile:
			line = line.strip().split('\t')
			contig, bait_start, bait_end, bait_orient = line[0], int(line[3])-1, int(line[4]), line[5]
			bait_id = line[6].split('ID=')[-1].split(';')[0]
			target_start, target_end, target_orient = int(line[10])-1, int(line[11]), line[13]
			target_id = line[-1].split('ID=')[-1].split(';')[0]

			if bait_id == target_id and bait_start == target_start and bait_end == target_end:
				continue

			dist = -sys.maxsize
			if target_start < bait_end and bait_start < target_end:
				dist = 0
			elif target_orient == '+' and bait_start >= target_end:
				dist = bait_start - target_end + 1
			elif target_orient == '+' and bait_end <= target_start:
				dist = (target_start - bait_end + 1)*-1
			elif target_orient == '-' and bait_start >= target_end:
				dist = (bait_start - target_end + 1)*-1
			elif target_orient == '-' and bait_end <= target_start:
				dist = target_start - bait_end + 1
			elif bait_start >= target_end:
				dist = bait_start - target_end + 1
			elif bait_end <= target_start:
				dist = target_start - bait_end + 1
				
			clen = contig_length(contig)

			rpt_seq = crispr_info[bait_id]['rpt_seq']
			rpt_seq_len = len(rpt_seq)
			spacers = crispr_info[bait_id]['spacers']
			num_spacers = len(spacers)
			num_uniq_spacers = len(set(spacers))
			spacer_lengths = list(map(len, spacers))
			mean_spacer_length = np.mean(spacer_lengths)
			min_spacer_length = np.min(spacer_lengths)
			max_spacer_length = np.max(spacer_lengths)
			sdev_spacer_length = np.std(spacer_lengths)

			out.append([target_id[:24], bait_id, contig, clen, target_start, target_end, target_orient, bait_start, bait_end, bait_orient, dist, 'CRISPR', rpt_seq_len, rpt_seq, num_spacers, num_uniq_spacers, mean_spacer_length, min_spacer_length, max_spacer_length, sdev_spacer_length, ';'.join(spacers)])

	cur.close()
	conn.close()
	df = pd.DataFrame(out)
	colnames = ['target_id', 'bait_id', 'contig', 'contig_length', 'target_start', 'target_end', 'target_orient', 'bait_start', 'bait_end', 'bait_orient', 'dist', 'bait_type', 'rpt_seq_len', 'rpt_seq', 'num_spacers', 'num_uniq_spacers', 'mean_spacer_length', 'min_spacer_length', 'max_spacer_length', 'sdev_spacer_length', 'spacer_seqs']
	try:
		df.columns = colnames
	except:
		df = pd.DataFrame(columns=colnames)
	return df


def crispr_neighbors(input, output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('crispr_neighbors', conn)
		conn.close()

		print(db_path)
	
		gbucket_crispr_neighbors = join("crispestdb/" + db_path, basename(output.crispr_neighbors))
		if gbucket_exists(gbucket_crispr_neighbors) and compute_log_done and params.rerun is False:
			print("Downloading existing CRISPR neighbors...")
			gbucket_download(gbucket_crispr_neighbors, output.crispr_neighbors)
			return
	
		print("Running CRISPR Neighbors...")
		outdir = output.crispr_neighbors + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)

		print("Getting cas genes...")
		cas_genes = dict()
		with gzip.open(input.hmmsearch_cct, 'rt') as infile:
			infile.readline()
			for line in infile:
				line = line.strip().split('\t')
				evalue = float(line[2])
				score = float(line[3])
				if float(line[2]) < 1e-6 and '|Cas' in line[1]:
					if line[0] not in cas_genes:
						cas_genes[line[0]] = (line[1], evalue, score)
					else:
						if evalue < cas_genes[line[0]][1]:
							cas_genes[line[0]] = (line[1], evalue, score)
				
		print("Sorting prodigal GFF...")
		tmp_genes_gff = join(outdir, 'all_genes.gff')
		shell("bedtools sort -i {gff} > {out}".format(gff=input.gff, out=tmp_genes_gff))

		print("Counting Cas genes...")
		tmp_gff = join(outdir, 'cas_genes.gff')
		with open(tmp_gff, 'w') as outfile:
			with open(tmp_genes_gff) as infile:
				for line in infile:
		
					if line.startswith('#'):
						continue
		
					prot = line.strip().split('\t')[-1].split(';')[0].replace('ID=', '')
					if prot in cas_genes:
						outfile.write(line)

		tmp_cas_intersect = join(outdir, 'cas_intersect.gff')
		cmd = """cat {cas_genes} | awk -F'\t' '{{print $1"\t"$4-25000"\t"$5+25000"\t"$4"\t"$5"\t"$7"\t"$9}}' | awk -F'\t' '{{if ($2 < 0){{print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}else{{print $0}}}}' | bedtools intersect -sorted -a stdin -b {all_genes} -wa -wb > {out}""".format(cas_genes=tmp_gff, all_genes=tmp_genes_gff, out=tmp_cas_intersect)
		print(cmd)
		os.system(cmd)

		cas_neighbors = make_cas_neighbors_df(tmp_cas_intersect, cas_genes, input.sqldb)

		print("Writing minced gff file...")
		crispr_info = defaultdict(lambda : {'rpt_seq':None, 'spacers':[]})
		tmp_minced = join(outdir, 'minced.gff')
		outfile = open(tmp_minced, 'w')
		with gzip.open(input.minced, 'rt') as infile:
			for line in infile:
				line = line.strip().split()
				if line[2] == 'spacer':
					crispr = line[-1].split('Parent=')[-1].split(';')[0]
					crispr_info[crispr]['spacers'].append(line[-1].split('seq=')[-1].split(';')[0])
				if line[2] != 'repeat_region':
					continue
				crispr = line[-1].split('ID=')[-1].split(';')[0]
				rpt_seq = line[-1].split('rpt_unit_seq=')[-1].split(';')[0]
				crispr_info[crispr]['rpt_seq'] = rpt_seq
				print(*line, sep='\t', file=outfile)
		outfile.close()

		tmp_minced_sorted = join(outdir, 'minced.sorted.gff')
		shell("bedtools sort -i {minced} > {minced_srt}".format(minced=tmp_minced, gff=tmp_genes_gff, minced_srt=tmp_minced_sorted))

		tmp_minced_intersect = join(outdir, 'minced_intersect.gff')
		cmd = """cat {minced_srt} | awk -F'\t' '{{print $1"\t"$4-25000"\t"$5+25000"\t"$4"\t"$5"\t"$7"\t"$9}}' | awk -F'\t' '{{if ($2 < 0){{print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}}else{{print $0}}}}' | bedtools intersect -sorted -a stdin -b {all_genes} -wa -wb > {out}""".format(minced_srt=tmp_minced_sorted, all_genes=tmp_genes_gff, out=tmp_minced_intersect)
		os.system(cmd)

		crispr_neighbors = make_crispr_neighbors_df(tmp_minced_intersect, crispr_info, input.sqldb)

		neighbor_columns = ['target_id', 'bait_id', 'contig', 'contig_length', 'target_start', 'target_end', 'target_orient', 'bait_start', 'bait_end', 'bait_orient', 'dist', 'bait_type']
		cas_neighbors_only = cas_neighbors[neighbor_columns].drop_duplicates()
		crispr_neighbors_only = crispr_neighbors[neighbor_columns].drop_duplicates()

		neighbors_only = pd.concat([cas_neighbors_only, crispr_neighbors_only])
		cas_bait_only = cas_neighbors[['bait_id', 'bait_domain', 'bait_evalue', 'bait_score']].drop_duplicates()
		crispr_bait_only = crispr_neighbors[['bait_id', 'rpt_seq_len', 'rpt_seq', 'num_spacers', 'num_uniq_spacers', 'mean_spacer_length', 'min_spacer_length', 'max_spacer_length', 'sdev_spacer_length', 'spacer_seqs']].drop_duplicates()

		tmp_neighbors = join(outdir, 'neighbors.csv')
		tmp_cas_bait = join(outdir, 'cas.csv')
		tmp_crispr_bait = join(outdir, 'crispr.csv')

		conn = sqlite3.connect(output.crispr_neighbors)
		cur = conn.cursor()
		cur.executescript("""
			CREATE TABLE neighbors (
				target_id TEXT,
				bait_id TEXT,
				contig TEXT,
			 	contig_length INTEGER, 
				target_start INTEGER, 
				target_end INTEGER,
				target_orient TEXT, 
				bait_start INTEGER,
				bait_end INTEGER, 
				bait_orient TEXT,
				dist INTEGER,
				bait_type TEXT	
			);
			CREATE INDEX neighbors_target_id_idx ON neighbors(target_id);
			CREATE INDEX neighbors_bait_id_idx ON neighbors(bait_id);
			CREATE INDEX neighbors_contig_idx ON neighbors(contig);

			CREATE TABLE cas (
				bait_id TEXT,
				bait_domain TEXT,
				bait_evalue NUMERIC,
				bait_score NUMERIC	
			);
			CREATE INDEX cas_bait_id_idx ON cas(bait_id);
			CREATE TABLE crispr (
				bait_id TEXT,
				rpt_seq_len INTEGER,
				rpt_seq TEXT,
				num_spacers INTEGER,
				num_uniq_spacers INTEGER,
				mean_spacer_length NUMERIC,
				min_spacer_length INTEGER,
				max_spacer_length INTEGER,
				sdev_spacer_length NUMERIC,
				spacer_seqs TEXT
			);
			CREATE INDEX crispr_bait_id_idx ON crispr(bait_id);
		""")
		
		neighbors_only.to_sql('neighbors', conn, if_exists='append', index=False)
		cas_bait_only.to_sql('cas', conn, if_exists='append', index=False)
		crispr_bait_only.to_sql('crispr', conn, if_exists='append', index=False)

		cur.close()
		conn.close()

		gbucket_upload(output.crispr_neighbors, gbucket_crispr_neighbors)

		conn = get_connection()
		compute_log_complete('crispr_neighbors', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def make_sqldb(output, params):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		print(db_path)
		compute_log_done = get_compute_log('make_sqldb', conn)
		conn.close()
	
		gbucket_sqldb = join("crispestdb/" + db_path, basename(output.sqldb))
		if gbucket_exists(gbucket_sqldb) and compute_log_done and params.rerun is False:
			print("Downloading existing SQL database...")
			gbucket_download(gbucket_sqldb, output.sqldb)
			return None

		print("Making SQL database...")
		outdir = output.sqldb + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)
		
		gbucket_neighbors = join("crispestdb/" + db_path, params.sample+'.crispr_neighbors.tsv.gz')
		gbucket_gff = join("crispestdb/" + db_path, params.sample+'.prodigal.gff.gz')
		gbucket_contigs = join("crispestdb/" + db_path, params.sample+'.contigs.tsv.gz')

		out_neighbors = join(outdir, basename(gbucket_neighbors))
		out_gff = join(outdir, basename(gbucket_gff))
		out_contigs = join(outdir, basename(gbucket_contigs))

		gbucket_download(gbucket_neighbors, out_neighbors)
		gbucket_download(gbucket_gff, out_gff)
		gbucket_download(gbucket_contigs, out_contigs)

		conn = sqlite3.connect(output.sqldb)
		cur = conn.cursor()
		cur.executescript("""
			CREATE TABLE gene_coords(
				contig TEXT NOT NULL,
				start INTEGER NOT NULL,
				end INTEGER NOT NULL,
				orient TEXT NOT NULL,
				protein_id TEXT NOT NULL
			);
			CREATE TABLE contig_lengths(
				contig TEXT NOT NULL,
				length INTEGER NOT NULL
			);
			CREATE TABLE crispr_neighbors(
				contig TEXT NOT NULL,
				protein_id TEXT NOT NULL,
				dist INTEGER NOT NULL,
				info TEXT NOT NULL
			);
			CREATE INDEX gene_coords_contig_idx ON gene_coords(contig);
			CREATE INDEX gene_coords_protein_id_idx ON gene_coords(protein_id);
			CREATE INDEX crispr_neighbors_contig_idx ON crispr_neighbors(contig);
			CREATE INDEX crispr_neighbors_protein_id_idx ON crispr_neighbors(protein_id);
			CREATE UNIQUE INDEX contig_lengths_contig_idx ON contig_lengths(contig);
		""")

		with gzip.open(out_neighbors, 'rt') as infile:
			infile.readline()
			for line in infile:
				line = line.strip().split('\t')
				cur.execute('INSERT INTO crispr_neighbors(contig,protein_id,dist,info) VALUES(?,?,?,?)', line)

		with gzip.open(out_gff, 'rt') as infile:
			for line in infile:
				if line.startswith("#"):
					continue
				line = line.strip().split('\t')
				contig, start, end, orient, pid = line[0], line[3], line[4], line[5], line[-1].replace('ID=', '').split(';')[0]
				line = [contig, start, end, orient, pid]
				cur.execute('INSERT INTO gene_coords(contig,start,end,orient,protein_id) VALUES(?,?,?,?,?)', line)

		with gzip.open(out_contigs, 'rt') as infile:
			infile.readline()
			for line in infile:
				line = line.strip().split('\t')	
				cur.execute('INSERT INTO contig_lengths(contig,length) VALUES(?,?)', [line[0], line[-1]])

		conn.commit()
		conn.close()

		gbucket_upload(output.sqldb, gbucket_sqldb)

		conn = get_connection()
		compute_log_complete('make_sqldb', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def blast_casdelta_flanks(input, output, params):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('blast_casdelta_flanks', conn)
		conn.close()

		gbucket_blast_info = join("crispestdb/" + db_path, basename(db_path), basename(output.blast_info))
		gbucket_regions_fna = join("crispestdb/" + db_path, basename(db_path), basename(output.regions_fna))
		if gbucket_exists(gbucket_blast_info) and gbucket_exists(gbucket_regions_fna) and compute_log_done and params.rerun is False:
			print("Downloading existing casdelta_flanks...")
			gbucket_download(gbucket_blast_info, output.blast_info)
			gbucket_download(gbucket_regions_fna, output.regions_fna)
			return

		outdir = output.blast_info + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)

		gbucket_flanks = join("crispestdb-other/casdelta_flanks.fna")
		flanks_out = join(outdir, basename(gbucket_flanks))
		gbucket_download(gbucket_flanks, flanks_out)
		
		blast_out = join(outdir, 'flank_blast.tsv')
		cmd = "blastn -task blastn -evalue 1e-12 -query {query_fna} -db {blast_pref} -outfmt '6 qseqid sseqid length qlen slen qstart qend sstart send pident evalue bitscore qseq sseq' 1> {out} 2> /dev/null".format(query_fna=flanks_out, out=blast_out, blast_pref=input.blastdb)
		print("command:", cmd)
		shell(cmd)

		flanks_columns = ['qseqid', 'sseqid', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'pident', 'evalue', 'bitscore', 'qseq', 'sseq']
		try:
			flanks = pd.read_csv(blast_out, sep='\t', header=None)
			flanks.columns = flanks_columns
		except:
			flanks = pd.DataFrame(columns=flanks_columns)

		flanks = flanks[(flanks.qend - flanks.qstart) / flanks.qlen > 0.95]
		flanks = flanks[flanks.pident > 95.0]

		flanks['p100_id'] = [q.split('_')[0] for q in flanks.qseqid]
		flanks['end'] = [q.split('_')[1] for q in flanks.qseqid]
		flanks['num'] = [q.split('_')[2] for q in flanks.qseqid]

		flanks['sseqid'] = [s.split('|')[1] if s.count('|') >= 2 else s for s in flanks.sseqid]
		unique_flank_pairs = set([(pid, num) for pid, num in zip(flanks.p100_id, flanks.num)])

		genome_fasta = pyfastx.Fasta(input.genome)

		seq_fragments = dict()
		seq_fragment_info = dict()
		for pid, num in unique_flank_pairs:

			flanks_filt = flanks[(flanks.p100_id==pid) & (flanks.num == num)]
			flanks_filt_left = flanks_filt[flanks_filt.end == 'left']
			flanks_filt_right = flanks_filt[flanks_filt.end == 'right']

			if flanks_filt_left.shape[0] == 0 or flanks_filt_right.shape[0] == 0:
				continue

			for contig, sstart, send in zip(flanks_filt_left.sseqid, flanks_filt_left.sstart, flanks_filt_left.send):

				flanks_filt_left_match = flanks_filt_left[(flanks_filt_left.sseqid == contig) & (flanks_filt_left.sstart == sstart) & (flanks_filt_left.send == send)]
				flanks_filt_right_match = flanks_filt_right[flanks_filt_right.sseqid == contig]

				for rcontig, rstart, rend in zip(flanks_filt_right_match.sseqid, flanks_filt_right_match.sstart, flanks_filt_right_match.send):
					flanks_filt_right_match_single = flanks_filt_right_match[(flanks_filt_right_match.sseqid == rcontig) & (flanks_filt_right_match.sstart == rstart) & (flanks_filt_right_match.send == rend)]
					
					left_orient, right_orient = 'F', 'F'
					if send < sstart:
						left_orient = 'R'
					if rend < rstart:
						right_orient = 'R'
				
					if left_orient != right_orient:
						continue

					full_start, full_end = min([send, sstart, rend, rstart]), max([send, sstart, rend, rstart])
					length = full_end - full_start

					full_seq = str(genome_fasta.fetch(contig, (full_start, full_end))).upper()

					if left_orient == 'R':
						full_seq = revcomp(full_seq)

					name = contig +'__'+str(full_start)+'__'+str(full_end)+'__'+pid+'__'+num
					seq_fragments[name] = full_seq

					seq_fragment_info[name] = (list(flanks_filt_left_match.length)[0], list(flanks_filt_left_match.pident)[0], list(flanks_filt_right_match_single.length)[0], list(flanks_filt_right_match_single.pident)[0], left_orient, contig, full_start, full_end, len(genome_fasta[contig]))

		outfile = open(output.regions_fna, 'w')
		for n in seq_fragments:
			print('>'+n, seq_fragments[n], sep='\n', file=outfile)
		outfile.close()

		outfile = open(output.blast_info, 'w')
		print("seqid", "left_length", "left_pident", "right_length", "right_pident", "orient", "contig", "start", "end", "contig_length", sep='\t', file=outfile)
		for n in seq_fragment_info:
			print(n, *seq_fragment_info[n], sep='\t', file=outfile)
		outfile.close()

		gbucket_upload(output.blast_info, gbucket_blast_info)
		gbucket_upload(output.regions_fna, gbucket_regions_fna)

		conn = get_connection()
		compute_log_complete('blast_casdelta_flanks', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


def annotate_casdelta_flanks(input, output, params):
	
	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('annotate_casdelta_flanks', conn)
		conn.close()
	
		gbucket_gff = join("crispestdb/" + db_path, basename(output.gff))
		if gbucket_exists(gbucket_gff) and compute_log_done and params.rerun is False:
			print("Downloading existing annotation files...")
			gbucket_download(gbucket_gff, output.gff)
			return
		
		outdir = join(pwd, output.gff + '.outdir')
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)
	
		tmp_fna = join(outdir, 'tmp.fna')
		out_tmp_fna = open(tmp_fna, 'w')
		contig_map = dict()
		for i, rec in enumerate(SeqIO.parse(input.regions_fna, "fasta")):
			contig_map[str(i+1)] = rec.id
			rec.id = str(i+1)
			SeqIO.write([rec], out_tmp_fna, 'fasta')
		out_tmp_fna.close()	

		if len(contig_map) == 0:
			
			o = open(output.gff, 'w')
			o.close()

			gbucket_upload(output.gff, gbucket_gff)
			conn = get_connection()
			compute_log_complete('annotate_casdelta_flanks', params, conn)
			conn.commit()
			return

		prokka_dir = join(outdir, 'prokka')
		prokka_prefix = 'prokka'

		cmd = 'singularity exec --bind /var:/var /home/mdurrant/prokka.sif prokka --metagenome --outdir {o} --prefix {p} {f}'.format(f=tmp_fna, o=prokka_dir, p=prokka_prefix)
		print('command:', cmd)
		shell(cmd)

		prokka_gff = join(prokka_dir, 'prokka.gff')
		prokka_faa = join(prokka_dir, 'prokka.faa')
		id2pid = dict()
		for rec in SeqIO.parse(prokka_faa, 'fasta'):
			id2pid[rec.id] = seqhash(str(rec.seq), seqtype='prot')

		gff_recs = []
		with open(prokka_gff) as infile:
			for line in infile:
				if line.startswith("#"):
					continue
				line = line.strip().split('\t')
				if len(line) != 9:
					continue
				line[0] = contig_map[line[0]]

				if line[2] == "CDS":
					sid = line[-1].split('ID=')[-1].split(';')[0]
					pid = id2pid[sid]
					line[-1] = re.sub(sid, pid, line[-1])
				gff_recs.append(line)
				
		outminced = run_minced_simple(tmp_fna, outdir)
		with open(outminced) as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				line = line.strip().split('\t')
				if len(line) != 9:
					continue
				line[0] = contig_map[line[0]]
				gff_recs.append(line)

		outccf = run_crisprcasfinder_simple(tmp_fna, outdir, add_pwd=False)
		with open(outccf) as infile:
			for line in infile:
				if line.startswith('#'):
					continue
				line = line.strip().split('\t')
				if len(line) != 9:
					continue
				line[0] = contig_map[line[0]]
				gff_recs.append(line)

		outfile = open(output.gff, 'w')
		for r in gff_recs:
			print(*r, sep='\t', file=outfile)
		outfile.close()
		
		gbucket_upload(output.gff, gbucket_gff)
		conn = get_connection()
		compute_log_complete('annotate_casdelta_flanks', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)
	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()

def cmsearch_omegarna(input, output, params, threads):

	conn = None
	outdir = None
	try:
		conn = get_connection()
		db_path = get_sample_field('db_path', params, conn)
		compute_log_done = get_compute_log('cmsearch_omegarna', conn)
		conn.close()	

		gbucket_cmsearch_out = join("crispestdb/" + db_path, basename(output.cmsearch_out))

		if gbucket_exists(gbucket_cmsearch_out) and compute_log_done and params.rerun is False:
			print("Downloading existing cmsearch omegarna files...")
			gbucket_download(gbucket_cmsearch_out, output.cmsearch_out)
			return

		outdir = output.cmsearch_out + '.outdir'
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		makedirs(outdir, exist_ok=True)
		#makedirs('data', exist_ok=True)
		#if not isfile('data/Koonin_IscB_IsrB_omegaRNA.cm'):
		#	gbucket_omegarna_profiles = join("crispestdb-other", 'Koonin_IscB_IsrB_omegaRNA.cm')
		#	gbucket_download(gbucket_omegarna_profiles, 'data/Koonin_IscB_IsrB_omegaRNA.cm')

		print("Gunzipping genome...")
		gunzipped_genome = gunzip_file(input.genome, outdir)
		outtbl = join(outdir, basename(output.cmsearch_out.replace('.tsv', '')+'.tbl'))

		print("Running cmsearch...")
		command = 'cmsearch -E 10.0 --noali --cpu {threads} --tblout {tbl} /home/mdurrant/data/Koonin_IscB_IsrB_omegaRNA.cm {genome}'.format(tbl=outtbl, genome=gunzipped_genome, threads=threads)
		print('cmsearch command:', command)
		shell(command)

		outfile = open(output.cmsearch_out, 'w')
		header = ['target_name', 'target_accession', 'query_name', 'query_accession', 'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'evalue', 'inc', 'desc']
		print(*header, sep='\t', file=outfile)
		with open(outtbl) as infile:
			for line in infile:
				if line.startswith('#'): continue
				line = line.strip().split()
				
				print(*line[:len(header)-1], ' '.join(line[len(header)-1:]), sep='\t', file=outfile)
		outfile.close()

		gbucket_upload(output.cmsearch_out, gbucket_cmsearch_out)
		conn = get_connection()
		compute_log_complete('cmsearch_omegarna', params, conn)
		conn.commit()

	except Exception as e:
		raise(e)

	finally:
		if outdir is not None and os.path.isdir(outdir):
			shutil.rmtree(outdir)
		if conn is not None:
			conn.close()


if __name__ == '__main__':
	
	if sys.argv[1] == 'download_genome':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		output = OutputFiles(fromdict={'genome': args['genome']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.genome), exist_ok=True)
		download_genome(output, params)

	if sys.argv[1] == 'pyfastx_genome':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome']})
		output = OutputFiles(fromdict={'genome_pyfastx':args['genome_pyfastx']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.genome_pyfastx), exist_ok=True)
		pyfastx_genome(input, output)

	if sys.argv[1] == 'genome_stats':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome']})
		output = OutputFiles(fromdict={'genome_stats': args['genome_stats'], 'contigs': args['contigs']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.genome_stats), exist_ok=True)
		genome_stats(input, output, params)

	if sys.argv[1] == 'prodigal':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome']})
		output = OutputFiles(fromdict={'faa': args['faa'], 'gff': args['gff']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.faa), exist_ok=True)
		prodigal(input, output, params, threads)
	
	if sys.argv[1] == 'hmmsearch_crisprcastyper':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'faa': args['faa'], 'gff': args['gff'], 'genome_stats': args['genome_stats']})
		output = OutputFiles(fromdict={'hmmsearch_out': args['hmmsearch_out']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.hmmsearch_out), exist_ok=True)
		hmmsearch_crisprcastyper(input, output, params, threads)
	
	if sys.argv[1] == 'mash_sketch':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome']})
		output = OutputFiles(fromdict={'mash_out': args['mash_out']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.mash_out), exist_ok=True)
		mash_sketch(input, output, params, threads)

	if sys.argv[1] == 'genomesearch_markers':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'faa': args['faa'], 'genome_stats': args['genome_stats']})
		output = OutputFiles(fromdict={'genomesearch_markers': args['genomesearch_markers']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.genomesearch_markers), exist_ok=True)
		genomesearch_markers(input, output, params, threads)
		
	if sys.argv[1] == 'minced':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome'], 'genome_stats':args['genome_stats']})
		output = OutputFiles(fromdict={'minced': args['minced']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.minced), exist_ok=True)
		minced(input, output, params, threads)

	if sys.argv[1] == 'crisprcasfinder':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome'], 'genome_stats':args['genome_stats']})
		output = OutputFiles(fromdict={'crisprcasfinder': args['crisprcasfinder']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.crisprcasfinder), exist_ok=True)
		crisprcasfinder(input, output, params, threads)

	if sys.argv[1] == 'find_cas_fusions':	
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'hmmsearch_out': args['hmmsearch_out'], 'faa': args['faa']})
		output = OutputFiles(fromdict={'cas_fusions_bed': args['cas_fusions_bed'], 'cas_fusions_summary':args['cas_fusions_summary'], 'crisprcastyper_domains':args['crisprcastyper_domains']})
		params = Params(fromdict={'crisprcastyper_hmm':args['crisprcastyper_hmm'], 'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.cas_fusions_bed), exist_ok=True)
		find_cas_fusions(input, output, params)

	if sys.argv[1] == 'hmmsearch_isescan':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'faa': args['faa']})
		output = OutputFiles(fromdict={'hmmsearch_isescan_out': args['hmmsearch_isescan_out']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True', 'isescan_hmm_path': args['isescan_hmm_path']})
		threads = int(args['threads'])
		makedirs(dirname(output.hmmsearch_isescan_out), exist_ok=True)
		hmmsearch_isescan(input, output, params, threads)

	if sys.argv[1] == 'makeblastdb':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome']})
		output = OutputFiles(fromdict={'blastdb_nhr': args['blastdb_nhr'], 'blastdb_nin': args['blastdb_nin'],'blastdb_nog': args['blastdb_nog'],'blastdb_nsd': args['blastdb_nsd'],'blastdb_nsi': args['blastdb_nsi'],'blastdb_nsq': args['blastdb_nsq']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.blastdb_nhr), exist_ok=True)
		makeblastdb(input, output, params)

	if sys.argv[1] == 'self_targeting_spacers':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome'], 'crisprcasfinder': args['crisprcasfinder'], 'minced': args['minced'], 'faa': args['faa'], 'gff': args['gff']})
		output = OutputFiles(fromdict={'self_targeting_spacers': args['self_targeting_spacers']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.self_targeting_spacers), exist_ok=True)
		self_targeting_spacers(input, output, params)

	if sys.argv[1] == 'make_sqldb':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		output = OutputFiles(fromdict={'sqldb': args['sqldb']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.sqldb), exist_ok=True)
		make_sqldb(output, params)

	if sys.argv[1] == 'crispr_neighbors':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome': args['genome'], 'hmmsearch_cct': args['hmmsearch_cct'], 'minced': args['minced'], 'faa': args['faa'], 'gff': args['gff'], 'sqldb': args['sqldb']})
		output = OutputFiles(fromdict={'crispr_neighbors': args['crispr_neighbors']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.crispr_neighbors), exist_ok=True)
		crispr_neighbors(input, output, params)

	if sys.argv[1] == 'blast_casdelta_flanks':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'blastdb': args['blastdb'], 'genome':args['genome']})
		output = OutputFiles(fromdict={'blast_info':args['blast_info'], 'regions_fna':args['regions_fna']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.blast_info), exist_ok=True)
		blast_casdelta_flanks(input, output, params)
		
	if sys.argv[1] == 'annotate_casdelta_flanks':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'regions_fna':args['regions_fna']})
		output = OutputFiles(fromdict={'gff':args['gff']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		makedirs(dirname(output.gff), exist_ok=True)
		annotate_casdelta_flanks(input, output, params)
		
	if sys.argv[1] == 'cmsearch_omegarna':
		args = {arg.split('=')[0]:arg.split('=')[1] for arg in sys.argv[2:]}
		input = InputFiles(fromdict={'genome':args['genome']})
		output = OutputFiles(fromdict={'cmsearch_out':args['cmsearch_out']})
		params = Params(fromdict={'sample': args['sample'], 'rerun': args['rerun'] == 'True'})
		threads = int(args['threads'])
		makedirs(dirname(output.cmsearch_out), exist_ok=True)
		cmsearch_omegarna(input, output, params, threads)
