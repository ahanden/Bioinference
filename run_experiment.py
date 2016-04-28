#!/bin/python

# Performs an experiment and offers DETAILED documentation

from time import strftime
import MySQLdb
import networkx as nx
import newScripts as ns
import sys
import yaml

# Returns the average degree of a graph
def avgDegree(G):
    return sum(len(G[n]) for n in G)/float(len(G))

# Reads gene symbols from a file
def readGeneFile(fname,conn):
    genes = {}
    cursor = conn.cursor()
    with open(fname,'r') as file:
        for line in file:
            symbol = line.strip()
            cursor.execute("SELECT entrez_id FROM genes.genes WHERE symbol = %s",[symbol])
            row = cursor.fetchone()
            if row:
                if row[0] in genes:
                    print "Warning: Found duplicate/conflicting symbols - %s and %s." % (symbol, genes[row[0]])
                else:
                    genes[row[0]] = symbol
            else:
                cursor.execute("SELECT entrez_id FROM genes.gene_synonyms WHERE symbol = %s",[symbol])
                row = cursor.fetchone()
                if row:
                    if row[0] in genes:
                        print "Warning: Found duplicate/conflicting symbols - %s and %s." % (symbol, genes[row[0]])
                    else:
                        genes[row[0]] = symbol
                else:
                    cursor.execute("SELECT entrez_id FROM genes.discontinued_genes WHERE discontinued_symbol = %s",[symbol])
                    row = cursor.fetchone()
                    if row:
                        if row[0] in genes:
                            print "Warning: Found duplicate/conflicting symbols - %s and %s." % (symbol, genes[row[0]])
                        else:
                            genes[row[0]] = symbol
                    else:
                        print "Warning: Unable to find an Entrez ID for gene symbol %s" % (symbol)

    return genes


# Read the experiment parameters
args = {}
with open(sys.argv[1],'r') as file:
    args = yaml.load(file)
print "----------------------------------------------"
print "Analysis performed on %s" % (strftime("%c"))
print

# Connect to the database (or read a file)
conn = None
try:
    conn = MySQLdb.connect(
        host='localhost',
        db=args['database']['name'],
        user=args['database']['user'],
        passwd=args['database']['password'])
    print "Used %s database for interactions." % (args['database']['name'])
except KeyError:
    sys.stderr.write("Error: You must give a database name, user, and password in the YAML file.\n")
    sys.exit(1)

# Read in the interactome
CI = nx.Graph()
if 'format' not in args['interactome']:
    sys.stderr.write("Error: You must give an interactome format (database, GML, TSV, or SIF).\n")
    sys.exit(1)
elif args['interactome']['format'] == 'database':
    min_pubs = args['interactome'].get('min_pubs', 0)
    CI = ns.readDB(conn, min_pubs=min_pubs)
else:
    if 'file' not in args['interactome']:
        sys.stderr.write("Error: You must provide an input file if you are not reading the CI from a database.\n")
        sys.exit(1)
    elif args['interactome']['format'] == 'GML':
        CI = ns.readGML(['interactome']['file'])
    elif args['interactome']['format'] == 'SIF':
        CI = ns.readSIF(args['interactome']['file'])
    elif args['interactome']['format'] == 'TSV':
        CI = ns.readTSV(args['interactome']['file'])
    else:
        sys.stderr.write("Error: %s is an unknown file format for reading the interactome - must be database, GML, SIF, or TSV\n" % (args['interactome']['format']))
        sys.exit(1)

CI = max(nx.connected_component_subgraphs(CI), key=len)
print "Interactome stats: %d nodes and %d edges using a publication threshold of %d" % (len(CI), len(CI.edges()), min_pubs)
print

# Conditionally write the interactome
if 'output'in args['interactome']:
    if 'file' not in args['interactome']['output']:
        sys.stderr.write("Error: You must define a file to save the interactome to.\n")
        sys.exit(1)
    elif 'format' not in args['interactome']['output']:
        sys.stderr.write("Error: You must define a file format for the interactome's output. (GML, SIF, or TSV)\n")
        sys.exit(1)
    elif args['interactome']['output']['format'] == 'GML':
        ns.writeGML(CI,args['interactome']['output']['file'])
    elif args['interactome']['output']['format'] == 'TSV':
        ns.writeTSV(CI,args['interactome']['output']['file'])
    elif args['interactome']['output']['format'] == 'SIF':
        ns.writeSIF(CI,args['interactome']['output']['file'])
    else:
        sys.stderr.write("Error: %s is an unknown file format for writing the interactome - must be GML, SIF, or TSV.\n" % (args['interactome']['output']['format']))
        sys.exit(1)

# Read genes of interest
gois = set(readGeneFile(args['gois']['file'],conn))
seeds = gois & set(CI)
print "%d unique genes of interest (GOIs) read from file. %d of these have interactions." % (len(gois), len(seeds))
print

# Output genes of interest
if 'label_file' in args['gois']:
    with open(args['gois']['label_file'],'w') as file:
        file.write("Gene\tisGOI\n")
        for gene in gois:
            file.write("%d\t1\n" % (gene))
if 'true_orphan_file' in args['gois']:
    with open(args['gois']['true_orphan_file'],'w') as file:
        file.write("GOIs with no interactions\n")
        if len(gois - seeds) == 0:
            file.write("(None - all GOIs had interactions)\n")
        for gene in gois - seeds:
            file.write("%d\n" % (gene))

# Expand the network
G = nx.Graph()
if 'method' not in args['expansion']:
    sys.stderr.write("Error: You must define an expansion method (pval or sp).\n")
    sys.exit(1)

elif args['expansion']['method'] == 'sp':
    prune = args['expansion'].get('prune', True)
    iters = args['expansion'].get('iters', 1000)

    print "Expanded network using shortest paths method"
    print "Expansion parameters: max_dist=%s, min_pubs=%d, filter=%s" % (max_dist, min_pubs, filter)

    G = ns.spGraph(seeds, CI, prune=prune, iters=iters)

elif args['expansion']['method'] == 'pval':
    max_p    = args['expansion'].get('max_p',    0.05)
    max_size = args['expansion'].get('max_size', None)
    verbose  = args['expansion'].get('verbose',  False)

    print "Expanded network using p-value method"
    print "Expansion parameters: max_p=%f, max_size=%d, verbose=%s" % (max_p, max_size, verbose)

    G = pvalGraph(goi,CI,max_p=0.05,max_size=None, verbose=False)

else:
    sys.stderr.write("Error: %s is an unrecognized expansion method - must be either pval or sp\n" % (args['expansion']['method']))
    sys.exit(1)

print "Expanded network stats: %d nodes and %d edges and average degree %f. %d/%d GOIs included in the network." % (len(G), len(G.edges()), avgDegree(G), len(gois & set(G)), len(gois))
print

# Write the expanded network to file
if 'output' in args['expansion']:
    if 'format' not in args['expansion']['output']:
        sys.stderr.write("Error: You must define a file format for saving the expanded network (GML, SIF, or TSV).\n")
        sys.exit(1)
    elif 'file' not in args['expansion']['output']:
        sys.stderr.write("Error: You must give a file name for saving the expanded network.\n")
        sys.exit(1)
    elif args['expansion']['output']['format'] == 'GML':
        ns.writeGML(G,args['expansion']['output']['file'])
    elif args['expansion']['output']['format'] == 'SIF':
        ns.writeSIF(G,args['expansion']['output']['file'])
    elif args['expansion']['output']['format'] == 'TSV':
        ns.writeTSV(G,args['expansion']['output']['file'])
    else:
        sys.stderr.write("Error: %s is not a recognized file format for saving the expanded network - must be SIF, TSV, or GML.\n" % (args['expansion']['output']['format']))
        sys.exit(1)

# Clustering
if 'clustering' in args:
    cd = {}

    if 'method' not in args['clustering']:
        sys.stderr.write("Error: You must define a clustering method (infomap or louvain).\n")
        sys.exit(1)
    elif args['clustering']['method'] == 'infomap':
        print "Clustering used the InfoMap method"
        cd = ns.infomapCluster(G)
    elif args['clustering']['method'] == 'louvain':
        print "Clustering used the Louvain method"
        cd = ns.louvainCluster(G)
    else:
        sys.stderr.write("Error: %s is not a recognized clustering method - must be either infomap or louvain.\n")
        sys.exit(1)

    print "%d total clusters found. The largest cluster has %d nodes, the average size is %d." % (
        len(cd),
        max([len(cd[x]) for x in cd]),
        sum([len(cd[x]) for x in cd]) / float(len(cd)))

# Writing clusters
if 'file' in args['clustering']:
    with open(args['clustering']['file'],'w') as file:
        file.write("Node\tCluster\n")
        for node in G.nodes(data=True):
            file.write("%d\t%d\n" % (node[0],node[1]['cluster']))

##################### CLEAN UP##################
with open('pvals.tsv','w') as file:
    file.write("Edge\tPvalue\n")
    for edge in G.edges(data=True):
        file.write("%d (INTERACTS) %d\t%f\n" % (edge[0],edge[1],edge[2]['pval']))
