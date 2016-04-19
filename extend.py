import sys
import newScripts
import networkx as nx
import MySQLdb

sys.stderr.write("Reading interactome...\n")
conn = MySQLdb.connect(
    host='localhost',
    db='interactions',
    user='adam',
    passwd='chanlab2016',
)

CI = newScripts.readDB(conn)

cursor = conn.cursor()

sys.stderr.write("Reading seed genes...\n")
genes = []
with open(sys.argv[1],'r') as file:
    for line in file:
        symbol = line.strip()
        cursor.execute("SELECT entrez_id FROM genes.genes WHERE symbol = %s",[symbol])
        row = cursor.fetchone()
        if row:
            genes.append(row[0])
        else:
            cursor.execute("SELECT entrez_id FROM genes.gene_synonyms WHERE symbol = %s",[symbol])
            row = cursor.fetchone()
            if row:
                genes.append(row[0])
            else:
                cursor.execute("SELECT entrez_id FROM genes.discontinued_genes WHERE discontinued_symbol = %s",[symbol])
                row = cursor.fetchone()
                if row:
                    genes.append(row[0])
                else:
                    sys.stderr.write("Warning: unable to find an Entrez ID for gene symbol %s\n" % symbol)

#print "\n".join([str(x) for x in genes])

G = newScripts.extendGraph(genes,CI)
print "Extended graph has %d genes and %d edges" % (len(G), len(G.edges()))
newScripts.writeTSV(G,sys.argv[2])
