exclude_directories = [".git", ".idea", "venv"]
sql_tables = ["proteins", "samples", "crisprs", "contigs", "contig2sample", "crisprcoords", "proteincoords", "prot2prot",
                 "prot2crispr", "prot2protwindow", "prot2crisprwindow"]
temp_files = ["contig2sample.tmp.sql.csv", "contigs.tmp.sql.csv", "crispr_coords.tmp.sql.csv", "CRISPRs.tmp.sql.csv",
                       "gene_coords.tmp.sql.csv", "merged_sorted_coords.tmp.sql.csv", "merged.sorted.tmp.sql.gff",
                       "protein2crispr_window.tmp.sql.csv", "protein2crispr.tmp.sql.csv", "protein2protein_window.tmp.sql.csv",
                       "protein2protein.tmp.sql.csv", "proteins.tmp.sql.csv", "temp.minced.sql.gff"]