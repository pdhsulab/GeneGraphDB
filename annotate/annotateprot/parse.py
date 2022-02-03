#!/usr/bin/env python
# coding: utf-8
# %%

# %%


import os
from Bio import SeqIO


# %%


infile_path = "out/hmmer.tab"
outfile_path = "OUTPUT/596adb7fa02d385b18.examples.gff"


# %%


def parse_cct_hmm(infile_path, outfile_path, targetid):
    with open(infile_path, "r") as infile, open(outfile_path, "a") as out_gff:
        next(infile)
        for line in infile:
            line_elem = line.split('\t')
            Hmm, ORF, tlen, qlen, Eval, score, start, end, Acc, Pos, Cov_seq, Cov_hmm, strand = line_elem
            strand, start, end = ".", 1, int(int(end)/ 3)
            seqname, source, feature, score, frame = ORF, "CCT", "CCT_domain", 1, 0
            seqname = targetid
            attrib = ';'.join(['ID=' + Hmm, 
                              'full_seq_evalue=' + str(Eval), 
                              'full_seq_score=' + str(score),
                              'qlen=' + str(qlen),
                              'tlen=' + str(tlen),
                              'qaln_perc=' + str(float(Cov_hmm)*100),
                              'taln_perc=' + str(float(Cov_seq)*100),
                             ])
            out = [seqname, source, feature, start, end, score, strand, frame, attrib]
     
            if float(Eval) < 1e-4:
                print(*out, sep = '\t', file=out_gff)



def parse_hmmsearch(infile_path, outfile_path):
    with open(infile_path, "r") as infile, open(outfile_path, "a") as out_gff:
        next(infile)
        for line in infile:
            line_elem = line.strip('\n').split(',')
            query_id,query_accession,query_len,hit_id,hit_score,hit_evalue,target_length,qstart,qend,tstart,tend = line_elem
            seqname, strand = hit_id, '.'
            source, feature, score, frame = "hmmsearch", "pfam_domain", "1", "0"
            attrib = ';'.join(['ID=' + query_id, 
                              'domain_evalue=' + str(hit_evalue), 
                              'domain_score=' + str(hit_score),
                              'qlen=' + str(query_len),
                              'tlen=' + str(target_length),
                               'qstart=' + str(qstart),
                              'qend=' + str(qend),
                              'tstart=' + str(tstart),
                               'tend=' + str(tend),
                             ])
            out = [seqname, source, feature, str(tstart), str(tend), score, strand, frame, attrib]
            if float(hit_evalue) < 1e-4:
                print(*out, sep = '\t', file=out_gff)


# %%


def put_gff_together(targetid):
    cct_out_path = 'cct/OUTPUT/' + targetid + '/hmmer.tab'
    gff_path = 'OUTPUT/' + targetid + '.gff'
    sorted_path = 'OUTPUT/' + targetid + '.sorted.gff'
    faa_path = 'INPUT/faa/' + targetid + '.faa'
    gff_path_final = 'OUTPUT/' + targetid + '.final.gff'
    hmm_outclean_path = 'hmmer/out_clean/' + targetid + '.csv'
    try:
        parse_cct_hmm(cct_out_path, gff_path, targetid)
    except:
        pass # some fna files impossible to retrieve at the moment
    parse_hmmsearch(hmm_outclean_path, gff_path)
    os.system("sortBed -i " + gff_path + " > " + sorted_path)
    os.system("cat " + sorted_path + " " + faa_path + " > " + gff_path_final)