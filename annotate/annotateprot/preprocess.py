from Bio import SeqIO
from Bio.Seq import Seq
import os
import sqlite3
import pandas as pd
import ast

infile_paths = ['INPUT/gff/' + file for file in os.listdir('INPUT/gff') if file != '.ipynb_checkpoints']
infile_paths

def get_prot_sequence(pid):
    con=sqlite3.connect("../../80kprotein_stats.db")
    cur = con.cursor()
    cmd = "SELECT sequence FROM proteins WHERE pid = '%s'" % pid 
    #print(cmd)
    cur.execute(cmd)
    return str(cur.fetchone()[0])
    con.close()

def get_permissive_rep(bait_pid):
    conn = sqlite3.connect('../../genegraph.db')
    cursor = conn.cursor()
    perm_rep = None
    cmd_p = "SELECT p30 FROM clusters WHERE p100 = '%s'" % (bait_pid)
    cursor.execute(cmd_p)
    try:
        perm_rep = cursor.fetchone()[0]
    except:
        pass
    #perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)

def create_gff(input_id):
    infile_path = "INPUT/gff/" + input_id + ".gff"
    sequence = get_prot_sequence(input_id)
    parent = get_permissive_rep(input_id)
    with open(infile_path, "w") as outfile:
        seqname, source, feature, start, end, score, strand, frame, attrib = input_id, 'PRODIGAL', 'CDS', '1', str(len(sequence)), '0', '.', '0', 'parent=' + parent
        header_gff = '\t'.join([seqname, source, feature, start, end, score, strand, frame, attrib])
        print(header_gff, file=outfile)
        print('>' + input_id, file=outfile)
        print(sequence, file=outfile)

# integrate later
def single_gff_to_faa(input_id):
    infile_path = "INPUT/gff/" + input_id + ".gff"
    multifasta_path = infile_path.replace("gff", "faa")
    start_reading = False
    outfile_gff_path = infile_path.replace("INPUT/gff", "OUTPUT")
    with open(infile_path, 'r') as infile, open(outfile_gff_path, 'w') as outfile_gff, open(multifasta_path, "w") as outfile:
        for line in infile:
            line = line.strip('\n')
            if line.startswith(">"):
                start_reading = True
            if not start_reading:
                print(line, file=outfile_gff)
            if start_reading:
                line = line.strip('*')
                print(line, file=outfile)

def create_fna(input_id):
    try:
        hit_near_tnpB_df = get_hit_near_inactive_tnpB(input_id, "../../tnpB_targetgenes_pfam.filtered.csv")
        hitid = list(hit_near_tnpB_df['hit_id'])[0]
        in_gff = "../ggdbfetch_output/" + hitid + ".examples.gff"
        out_fna = "INPUT/fna/" + input_id + ".fna"
        with open(in_gff, "r") as infile, open(out_fna, "w") as outfile:
            lines = infile.readlines()
            for line in lines:
                line = line.split('\t')
                start, end = int(line[3]), int(line[4])
                annot_id = line[-1].strip("ID=").strip("\n")
                if annot_id == input_id:
                    header = ">" + line[0]
                    break
            for i in range(len(lines)):
                line = lines[i].strip("\n")
                if line == header:
                    seq = lines[i+1][start-1:end]
                    print(header, file=outfile)
                    print(seq, file=outfile)
                    break
    except:
        pass

def get_hit_near_inactive_tnpB(inactive_tnpB, df_path):
    df = pd.read_csv(df_path).iloc[:,1:]
    baitlists_all = [ast.literal_eval(baitlist) for baitlist in df['baitp100s'].unique()]
    baitlists_expanded = []
    for baitlist in baitlists_all:
        for bait in baitlist:
            if bait == inactive_tnpB:
                baitlists_expanded.append(str(baitlist))
                break
    hits_with_inactive_tnpB_df = df[df["baitp100s"].isin(baitlists_expanded)].dropna()
    return hits_with_inactive_tnpB_df



# def single_gff_to_faa(input_id):
#     infile_path = "INPUT/gff/" + input_id + ".gff"
#     faa_path = infile_path.replace("gff", "faa")
#     fna_path = infile_path.replace("gff", "fna")
#     contigs = fasta_id_dict(fna_path)
#     #outfile_gff_path = infile_path.replace("INPUT/gff", "OUTPUT")
#     with open(infile_path, 'r') as infile, open(faa_path, "w") as outfile:
#         for line in infile:
#             line = line.strip('\n').split('\t')
#             if 'PRODIGAL' in line:
#                 header_old = line[0]
#                 header = ';'.join([line[0], line[3], line[4], line[6]])
#                 strand = line[6]
#                 try:
#                     dna_seq = contigs[header_old]
#                 except:
#                     pass #some proteins in contigs are annotated without reference sequence
#                 aa_seq = extract_faa_seq(header, dna_seq)
#                 print(aa_seq, header)
#                 asdf
#                 print('>' + header, file=outfile)
#                 print(aa_seq, file=outfile)
    

# def extract_faa_seq(header, seq):
#     header_list = header.split(';')
#     start = int(header_list[-3]) - 1
#     end = int(header_list[-2])
#     strand = header_list[-1]
#     if strand == '+':
#         dna = Seq(seq[start:end])
#         protein = dna.translate()
#         return protein
#     elif strand == '-':
#         dna = Seq(seq[start:end]).reverse_complement()
#         protein = dna.translate()
#         return protein
    
# def fasta_id_dict(fna_path):
#     contigs = {}
#     with open(fna_path) as handle:
#         for rec in SeqIO.parse(handle, 'fasta'):
#             header = rec.id
#             sequence = rec.seq
#             contigs[header] = sequence
#     return contigs