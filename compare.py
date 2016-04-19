import sys
import newScripts
import networkx as nx
import MySQLdb

conn = MySQLdb.connect(
    host='localhost',
    db='interactions',
    user='adam',
    passwd='chanlab2016',
)

old = newScripts.readGML(sys.argv[1])
new = newScripts.readDB(conn)

for node in new:
    new.node[node]['entrez'] = node

print "Old CI: %d nodes with %d edges" % (len(old), len(old.edges()))
print "New CI: %d nodes with %d edges" % (len(new), len(new.edges()))

old_nodes = set()
bad_nodes = set()
for node in old:
    try:
        old_nodes.update([int(old.node[node]['entrez'])])
    except:
        bad_nodes.update([old.node[node]['entrez']])

print "New CI contains %d of the old CI's genes" % (len(old_nodes & set(new)))
print "The old CI had %d bad node(s): %s" % (len(bad_nodes), ", ".join(bad_nodes))

count = 0
for edge in old.edges():
    try:
        if new.has_edge(int(old.node[edge[0]]['entrez']),int(old.node[edge[1]]['entrez'])):
            count += 1
    except:
        continue

print "The new CI contains %d of the old CI's edges" % (count)

cursor = conn.cursor()
count = 0
for node in old:
    try:
        cursor.execute("SELECT EXISTS(SELECT * FROM genes.genes WHERE entrez_id = %d)" % (int(old.node[node]['entrez'])))
        row = cursor.fetchone()
        if row[0]:
            count += 1
    except ValueError:
        continue
print "%d of the old CI's genes are still valid" % (count)
print old_nodes - set(new)
