#!/bin/python

# Performs an experiment and offers DETAILED documentation

from time import strftime
import MySQLdb
import networkx as nx
import signal
import sys
import yaml

from bioinference import graphIO, genes, expand, cluster

# Ctrl-C Interupt catching
def sigint_handler(signal, frame):
    print "SIGINT Detect: stopping analysis."
    sys.exit(1)

# Kill the program with an error
def die(message, status=1):
    sys.stderr.write("Error: %s\n" % message)
    sys.exit(status)

# Reads gene symbols from a file
def readGeneFile(fname,conn):
    genes = set()
    with open(fname,'r') as file:
        for line in file:
            symbol = line.strip()
            eids = genes.getEID(conn, symbol)
            if eids:
                genes.update(eids)
            else:
                print "Warning: Unable to find an Entrez ID for gene symbol %s" % (symbol)
    return genes

def run(arg):
    """Runs a pipeline of analyses based on a yaml configuration.
    See template.yaml for details.

    Parameters
    ----------
    yaml_file : path to the yaml configuration file

    OR

    parameters : a dictionary of paramters
    """

    signal.signal(signal.SIGINT, sigint_handler)

    # Format handlers
    writers = {
        'GML' : graphIO.writeGML,
        'TSV' : graphIO.writeTSV,
        'SIF' : graphIO.writeSIF
    }
    readers = {
        'GML' : graphIO.readGML,
        'SIF' : graphIO.readSIF,
        'TSV' : graphIO.readTSV
    }


    ############ READ YAML FILE ##############
    args = {}
    if type(arg) is dict:
        args = arg
    else:
        with open(arg,'r') as file:
            args = yaml.load(file)
    print "----------------------------------------------"
    print "Analysis performed on %s" % (strftime("%c"))
    print

    ############# CONNECT TO GENES DB ################
    genes_conn = None
    try:
        genes_conn = MySQLdb.connect(
            host='localhost',
            db=args['genes_db']['name'],
            user=args['genes_db']['user'],
            passwd=args['genes_db']['password'])
        print "Used %s database for genes." % (args['genes_db']['name'])
    except KeyError:
        die("You must give a database name, user, and password for the genes database.")

    ################ CONSOLIDATED INTERACTOME ################
    if 'interactome' not in args:
        die("You must provide an interactome as input")

    CI = nx.Graph()
    format = args['interactome'].get('format', None)
    if format is None:
        die("You must give an interactome format (database, GML, TSV, or SIF).")

    if format == 'database':
        min_pubs = args['interactome'].get('min_pubs', 0)

        try:
            conn = MySQLdb.connect(
                host='localhost',
                db=args['interactome']['ints_db']['name'],
                user=args['interactome']['ints_db']['user'],
                passwd=args['interactome']['ints_db']['password'])
            print "Used %s database for interactions." % (args['interactome']['ints_db']['name'])
            CI = graphIO.readDB(conn, min_pubs=min_pubs)
        except KeyError:
            die("You must give a database name, user, and password for the interactions database.")

    else:
        file  = args['interactome'].get('file', None)
        if file is None:
            die("You must provide an input file if you are not reading the CI from a database.")

        try:
            CI = readers[format](file)
        except KeyError:
            die("%s is an unknown file format for reading the interactome - must be database, GML, SIF, or TSV." % (format))

    CI = max(nx.connected_component_subgraphs(CI), key=len)
    print "Interactome stats: %d nodes and %d edges using a publication threshold of %d" % (len(CI), len(CI.edges()), min_pubs)
    print

    # Conditionally write the interactome
    if 'output'in args['interactome']:

        file = args['interactome']['output'].get('file',None)
        format = args['interactome']['output'].get('format',None)

        if file is None:
            die("You must define a file to save the interactome to.")
        if format is None:
            die("You must define a file format for the interactome's output - must be GML, SIF, or TSV.")

        try:
            writers[format](CI,file)
        except KeyError:
            die("%s is an unknown file format for writing the interactome - must be GML, SIF, or TSV." % (format))

    
    ############## GENES OF INTEREST ################
    gois = set()
    if 'gois' not in args:
        die("You must provide genes of interest as seeds.")
    if 'file' not in args['gois']:
        die("You must provide a file name containing genes of interest.")
    try:
        with open(args['gois']['file'],'r') as file:
            for line in file:
                symbol = line.strip()
                eids = genes.getEID(symbol, genes_conn)
                if eids:
                    gois.update(eids)
                else:
                    print "Warning: Unable to find an Entrez ID for gene symbol %s" % (symbol)
    except IOError as e:
        die("%s\nUnable to read GOI file." % (e))

    seeds = gois & set(CI)
    print "%d unique genes of interest (GOIs) read from file. %d of these have interactions." % (len(gois), len(seeds))
    print
    if len(seeds) == 0:
        print "Unable to continue if none of the GOIs have interactions."
        sys.exit(0)

    # Output genes of interest
    label_file  = args['gois'].get('label_file', None)
    orphan_file = args['gois'].get('true_orphan_file', None)
    if label_file is not None:
        with open(label_file,'w') as file:
            file.write("Gene\tisGOI\n")
            for gene in gois:
                file.write("%d\t1\n" % (gene))
    if orphan_file is not None:
        with open(orphan_file,'w') as file:
            file.write("GOIs with no interactions\n")
            if len(gois - seeds) == 0:
                file.write("(None - all GOIs had interactions)\n")
            for gene in gois - seeds:
                file.write("%d\n" % (gene))

    ########## NETWORK EXPANSION ###############
    if 'expansion' not in args:
        die("You must provide details on the expansion method to use.")
    G = nx.Graph()
    method = args['expansion'].get('method', None)
    if method is None:
        die("You must define an expansion method - must be sp or pval.")
    if method == 'sp':
        prune   = args['expansion'].get('prune', True)
        iters   = args['expansion'].get('iters', 1000)
        verbose = args['expansion'].get('verbose', False)
        max_p   = args['expansion'].get('max_p', 0.05)

        print "Expanded network using shortest paths method"
        print "Expansion parameters: prune=%s, iters=%s" % (prune, iters)

        if prune:
            G = expand.monteCarloPrune(seeds, CI, expand.linkerGraph, iters=iters, max_p=max_p)
        else:
            G = expand.linkerGraph(seeds, CI)

    elif method == 'pval':
        max_p    = args['expansion'].get('max_p',    0.05)
        max_size = args['expansion'].get('max_size', 2 * len(gois))
        verbose  = args['expansion'].get('verbose',  False)

        print "Expanded network using p-value method"
        print "Expansion parameters: max_p=%f, max_size=%d, verbose=%s" % (max_p, max_size, verbose)

        G = expand.cottrillGraph(gois,CI,max_p=0.05,max_size=None, verbose=verbose)

    else:
        die("%s is an unrecognized expansion method - must be either pval or sp." % (method))

    print "Expanded network stats: %d nodes and %d edges and average degree %f. %d/%d GOIs included in the network." % (
        len(G),
        len(G.edges()),
        sum(len(G[n]) for n in G)/float(len(G)),
        len(gois & set(G)),
        len(gois))
    print

    # Write the expanded network to file
    if 'output' in args['expansion']:
        format    = args['expansion']['output'].get('format', None)
        file      = args['expansion']['output'].get('graph_file', None)
        pval_file = args['expansion']['output'].get('pval_file', None)

        if format is None:
            die("You must define a file format for saving the expanded network - must be GML, SIF, or TSV.")
        if file is None:
            die("You must give a file name for saving the expanded network.")

        try:
            writers[format](G,file)
        except KeyError:
            die("%s is not a recognized file format for saving the expanded network - must be SIF, TSV, or GML.\n" % (format))

        if pval_file is not None and args['expansion']['method'] == 'sp':
            with open(pval_file,'w') as file:
                file.write("Edge\tP-value\n")
                for edge in G.edges(data=True):
                    file.write("%d (INTERACTS) %d\t%f\n" % (edge[0], edge[1], edge[2]['pval']))

    ################ CLUSTERING #################
    if 'clustering' in args:
        cd = {} # Cluster dictionary

        method       = args['clustering'].get('method', None)
        cluster_file = args['clustering'].get('file', None)
        weight       = args['clustering'].get('weight', False)

        if method is None:
            die("You must define a clustering method - must be infomap or louvain."
            )
        if method == 'infomap':
            print "Clustering used the InfoMap method. Using weights = %s" % (weight)
            if weight:
                cd = cluster.infomap(G,weight_feature='publications')
            else:
                cd = cluster.infomap(G)
        elif method == 'louvain':
            print "Clustering used the Louvain method. Using weights = %s" % (weight)
            if weight:
                cd = cluster.louvain(G,weight='publications')
            else:
                cd = cluster.louvain(G)
        else:
            die("%s is not a recognized clustering method - must be either infomap or louvain." % (method))

        print "%d total clusters found. The largest cluster has %d nodes, the average size is %d." % (
            len(cd),
            max([len(cd[x]) for x in cd]),
            sum([len(cd[x]) for x in cd]) / float(len(cd)))

        # Write clusters
        if cluster_file is not None:
            with open(cluster_file,'w') as file:
                file.write("Node\tCluster\n")
                for node in G.nodes(data=True):
                    file.write("%d\t%d\n" % (node[0],node[1]['cluster']))

if __name__ == "__main__":
    run(sys.argv[1])