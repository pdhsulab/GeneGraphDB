import sqlite3

conn = sqlite3.connect('genegraph.db')
cursor = conn.cursor()
def get_permissive_rep(stringent_rep):
    cmd_p = "SELECT reppid FROM permissive WHERE pid = '%s'" % (stringent_rep)
    cursor.execute(cmd_p)
    perm_rep = cursor.fetchone()[0]
    return(perm_rep)

print("loading stringent...")
p100p90 = dict()
with open("../clusters/stringent.csv") as infile:
    next(infile)
    for line in infile:
        line = line.strip().split(',')
        p100 = line[1][:18]
        p90 = line[0][:18]
        p100p90[p100] = p90

print("writing to file...")
p100_fail = []
with open("../clusters/OUTPUT/complete_clusters.tsv", "w") as outfile, open("../clusters/OUTPUT/p100_fail_list", "w") as err_outfile:
    print("p100", "p90", "p30", sep='\t', file=outfile)
    for p100 in p100p90:
        try:
            p90 = p100p90[p100]
            p30 = get_permissive_rep(p90)
            print(p100, p90, p30, sep='\t', file=outfile)
        except:
            print(p100, p100p90[p100], sep = ",", file=err_outfile)

conn.close()


