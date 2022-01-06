import os
from csv import reader

exclude_directories = [".git", ".idea", "venv"]
sql_tables = ["proteins", "samples", "crisprs", "contigs", "contig2sample", "crisprcoords", "proteincoords", "prot2prot",
                 "prot2crispr", "prot2protwindow", "prot2crisprwindow"]
temp_files_list = ["contig2sample.tmp.sql.csv", "contigs.tmp.sql.csv", "crispr_coords.tmp.sql.csv", "CRISPRs.tmp.sql.csv",
                       "gene_coords.tmp.sql.csv", "merged_sorted_coords.tmp.sql.csv", "merged.sorted.tmp.sql.gff",
                       "protein2crispr_window.tmp.sql.csv", "protein2crispr.tmp.sql.csv", "protein2protein_window.tmp.sql.csv",
                       "protein2protein.tmp.sql.csv", "proteins.tmp.sql.csv", "temp.minced.sql.gff"]
                   # "contig2sample.tmp.csv", "contigs.tmp.csv", "crispr_coords.tmp.csv", "CRISPRs.tmp.csv",
                   # "gene_coords.tmp.csv", "merged_sorted_coords.tmp.csv", "merged.sorted.tmp.gff",
                   # "protein2crispr_window.tmp.csv", "protein2crispr.tmp.csv", "protein2protein_window.tmp.csv",
                   # "protein2protein.tmp.csv", "proteins.tmp.csv", "temp.minced.gff"]
path_drep = "../drep_genomes/OUTPUT/rep_genomes/"
drep_samples = {}
# for directory in os.listdir(path_drep):
#     drep_samples[directory] = path_drep + "/" + directory
for directory in os.listdir(path_drep):
    path_2 = path_drep + directory
    if os.path.isdir(path_2):
        for directory_2 in os.listdir(path_2):
            path_3 = path_2 + "/" + directory_2
            if os.path.isdir(path_3):
                for directory_3 in os.listdir(path_3):
                    path_4 = path_3 + "/" + directory_3 + "/"
                    if os.path.isdir(path_4):
                        for directory_samp in os.listdir(path_4):
                            samp_dir = path_4 + directory_samp
                            if os.path.isdir(samp_dir):
                                drep_samples[directory_samp] = path_4
                
sampleids_error = []
with open("ggdb_multisql_errorlog.csv", "r") as f:
    infile = reader(f)
    for line in infile:
        sampleid = line[2].split(' : ')[0].split('/')[-1]
        sampleids_error.append(sampleid)
keys = set(sampleids_error).intersection(set(drep_samples.keys()))

drep_samples_error = {k:drep_samples[k] for k in keys}
#drep_samples_error = {"EARTH_3300025421_6": "../drep_genomes/OUTPUT/rep_genomes/A011/B026/C033/"}
