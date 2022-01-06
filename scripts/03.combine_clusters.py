import sqlite3

def get_permissive_rep(stringent_rep):
    conn = sqlite3.connect('../GeneGraphDB/genegraph.db')
    cursor = conn.cursor()
    cmd_p = "SELECT reppid FROM permissive WHERE pid = '%s'" % (stringent_rep)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    conn.close()
    return(perm_rep)

print("loading stringent...")
p100p90 = dict()
with open("stringent_stats.csv") as infile:
    next(infile)
    for line in infile:
        line = line.strip().split(',')
        p100p90[line[1]] = line[0]

# print("loading permissive...")
# p90p30 = dict()
# with open("permissive_stats.csv") as infile:
#     for line in infile:
#         line = line.strip().split(',')
#         p90p30[line[1]] = line[0]

print("writing to file...")
p100_fail = []
with open("OUTPUT/complete_clusters.tsv", "w") as outfile:
    print("p100", "p90", "p30", sep='\t', file=outfile)
    for p100 in p100p90:
        try:
            p90 = p100p90[p100]
            p30 = get_permissive_rep(p90)
            print(p100, p90, p30, sep='\t', file=outfile)
        except:
            p100_fail.append(p100)

with open("OUTPUT/p100_fail_list", "w") as outfile:
    for p100 in p100_fail:
        try:
            print(p100, p100p90[p100], sep = ",", file=outfile)
        except:
            print(p100, file=outfile)

# print("loading stringent...")
# p100p90 = dict()
# with open("stringent_stats.csv") as infile:
#     next(infile)
#     for line in infile:
#         line = line.strip().split(',')
#         p100p90[line[1]] = line[0]

# print("loading permissive...")
# p90p30 = dict()
# with open("permissive_stats.csv") as infile:
#     for line in infile:
#         line = line.strip().split(',')
#         p90p30[line[1]] = line[0]

# print("writing to file...")
# p90_fail = []
# with open("OUTPUT/complete_clusters.tmp.tsv", "w") as outfile:
#     print("p90", "p30", sep='\t', file=outfile)
#     for p90 in p90p30:
#         try:
#             p30 = p90p30[p90]
#             print(p90, p30, sep='\t', file=outfile)
#         except:
#             p90_fail.append(p90)

# with open("OUTPUT/p90_fail_list", "w") as outfile:
#     for p90 in p90_fail:
#         print(p90, file=outfile)
