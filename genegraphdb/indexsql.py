import sqlite3
def index_main():
    print("Creating database: %s" % DBNAME)
    con = sqlite3.connect('genegraph.db')
    cur = con.cursor()
    cur.execute('''CREATE TABLE proteins (hashid text, length real)''')
    cur.execute('''CREATE TABLE crisprs (hashid text)''')
    cur.execute('''CREATE TABLE contigs (hashid text, length real)''')
    cur.execute('''CREATE TABLE contig2sample (contighashid text, sampleid text)''')
    cur.execute('''CREATE TABLE crisprcoords (crisprhash text, contighash text, start real, end real)''')
    cur.execute('''CREATE TABLE proteincoords (phash text, contighash text, start real, end real, orientation text)''')
    cur.execute('''CREATE TABLE prot2prot (p1hash text, p2hash text)''')
    cur.execute('''CREATE TABLE prot2crispr (p1hash text, crisprhash text)''')
    cur.execute('''CREATE TABLE prot2protwindow (p1hash text, p2hash text)''')
    cur.execute('''CREATE TABLE prot2crisprwindow (p1hash text, crisprhash text)''')
    con.close()
    #os.system("sqlite3 genegraph.db")

    print("Indexes created")

def index_helper():