from fisher import pvalue
from math import log
from uuid import uuid4
import MySQLdb
import community
import networkx as nx
import os
import random
import re
import subprocess
import sys

####### GENE ID LOOKUP FUNCTIONS #########
def checkEID(cursor, eid):
    """Returns a valid Entrez ID given an uncertain Entrez ID.

    Entrez IDs may, from time to time, be discontinued for many
    reasons. This method interfaces with the gene database to 
    check that a given Entrez ID is still valid, or to find what
    the current valid Entrez ID is for an old, discontinued ID.

    Paramters
    ---------
    cursor : a MySQL cursor to the gene database

    eid : An Entrez ID

    Returns
    -------
    valid_eids : A set of valid Entrez IDs. The list will be
                 empty if no valid IDs were found.
    """

    valid_eids = set()
    cursor.execute("SELECT EXISTS(SELECT * FROM genes WHERE entrez_id = %(eid)s)", {'eid': eid})
    if cursor.fetchone()[0] == 1:
        valid_eids.add(eid)
    else:
        cursor.execute("SELECT entrez_id FROM discontinued_genes WHERe discontinued_id = %(eid)s", {'eid': eid})
        valid_eids.update([row[0] for row in cursor.fetchall()])
   
    return valid_eids

def getEID(cursor, symbol):
    """Returns a valid Entrez ID given a gene symbol.

    Paramters
    ---------
    cursor : a MySQL cursor to the gene database

    symbol : a gene symbol to convert

    Returns
    -------
    eids : a set of Entrez IDs - which will be empty if no valid IDs were found.
    """

    # Start with a direct query
    args = {"symbol": symbol}
    cursor.execute("SELECT entrez_id FROM genes WHERE symbol = %(symbol)s", args)
    eids = [row[0] for row in cursor.fetchall()]

    # If that didn't work, search for a discotninued symbol
    if not eids:
        cursor.execute("SELECT entrez_id FROM discontinued_genes WHERE discontinued_symbol = %(symbol)s", args)
        eids = [row[0] for row in cursor.fetchall()]

        # if THAT didn't work, search for synonyms
        if not eids:
            cursor.execute("SELECT entrez_id FROM gene_synonyms WHERE symbol = %(symbol)s", args)
            eids = [row[0] for row in cursor.fetchall()]

    # Return whatever we found (if anything)
    return set(eids)

def getSymbol(cursor, eid):
    """Returns a valid gene symbol given an Entrez ID.

    Paramters
    ---------
    cursor : a MySQL cursor to the gene database

    eid : An Entrez ID

    Returns
    -------
    symbols : a set of symbols - which will be empty if no valid symbols were found.
    """

    valid_eids = checkEID(cursor, eid)
    symbols = set()
    for gene in valid_eids:
        cursor.execute("SELECT symbol FROM genes WHERE entrez_id = %(eid)s AND symbol IS NOT NULL", {'eid': gene})
        symbols.update([row[0] for row in cursor.fetchall()])

    return symbols
   
def crossQuery(cursor, db, id):
    """Returns valid Entrez IDs given a foreign database and identifier.

    Paramters
    ---------
    cursor : a MySQL cursor to the gene database

    db : the name of the external database

    id : the foregin identifier

    Returns
    -------
    eids : a list of matching Entrez IDs - which will be empty if no matching ids were found.
    """
    
    cursor.execute("SELECT entrez_id FROM gene_xrefs WHERE Xref_db = %(db)s AND Xref_id = %(id)s", {'db': db, 'id': id})
    return [row[0] for row in cursor.fetchall()]

####### NETWORK LOADING FUNCTIONS #########
def readDB(conn,min_pubs=0):
    """Reads interactions into a NetworkX graph from a MySQL database

    Parameters
    ----------
    conn : A MySQLdb connection

    min_pubs : Minimum publications an edge must have in order to be included in the interactome (default=0)

    Returns
    -------
    G : A NetworkX Graph
    """
    G = nx.Graph()
    cursor = conn.cursor()
    cursor.execute("SELECT entrez_id1, entrez_id2 FROM interactions")
    cursor.execute("SELECT MAX(counted) FROM (SELECT COUNT(pubmed_id) AS counted FROM publications GROUP BY int_id) AS x")
    max_pubs = cursor.fetchone()[0]
    cursor.execute("SELECT entrez_id1, entrez_id2, count(pubmed_id) FROM interactions LEFT JOIN publications ON interactions.int_id = publications.int_id GROUP BY interactions.int_id HAVING count(pubmed_id) >= %s", [min_pubs])
    for row in cursor.fetchall():
        if(int(row[2]) >= min_pubs):
            G.add_edge(int(row[0]),int(row[1]),publications=int(row[2]),weight=1 + max_pubs - row[2])
    cursor.close()
    return G

def readGML(filename):
    """Reads a Graph from a GML formatted file.

    Parameters
    ----------
    filename : A file location

    Returns
    -------
    G : A NetworkX Graph
    """
    
    G = nx.read_gml(filename)
    # Remove self-edges
    for node in G.nodes():
        try:
            G.remove_edge(node,node)
        except:
            continue
    return G
    
def readTSV(filename):
    """Reads a Graph from a TSV formatted file.
    The TSV format should be as follows:
        EntrezID1    EntrezID2

    Parameters
    ----------
    filename : A file location

    Returns
    -------
    G : A NetworkX Graph
    """
    
    G = nx.Graph()
    with open(filename,'r') as file:
        header = True
        for line in file:
            if header:
                header = False
                continue
            (idA, idB) = line.strip().split("\t")
            if idA != idB: # Do not include self-interactions
                G.add_edge(idA, idB)
    return G
    
def readSIF(filename):
    """Reads a Graph from a SIF formatted file.

    Parameters
    ----------
    filename : A file location

    Returns
    -------
    G : A NetworkX Graph
    """
    
    G = nx.Graph()
    with open(filename,'r') as file:
        delimiter = " "
        
            
        for line in file:
            if "\t" in line:
                delimiter = "\t"

            fields = line.strip().split(delimiter)

            if len(fields) < 3:
                raise IndexError('The SIF file is improperly formatted. Visit http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats')
            
            eid1 = int(fields.pop(0))
            type = fields.pop(0)
            
            for eid2 in fields:
                if eid1 != eid2:
                    eid2 = int(eid2)
                    G.add_edge(eid1, eid2, type=type)
    return G

####### NETWORK WRITING FUNCTIONS ######
def writeGML(G, filename):
    """Writes a given network to file in GML format

    Parameters
    ----------
    G : A NetworkX graph
    
    filename : A file location
    """
    
    nx.write_gml(G,filename)

def writeTSV(G, filename):
    """Writes a given network to file in TSV format.
    The file's header will be "EntrezID1    EntrezID2"

    Parameters
    ----------
    G : A NetworkX graph
    
    filename : A file location
    """
    
    with open(filename, 'w') as file:
        file.write("EntrezID1\tEntrezId2\n")
        for edge in G.edges():
            file.write("%s\t%s\n" % edge)

def writeSIF(G, filename, compressed=False):
    """Writes a given network to file in SIF format.
    If an edge type is not found for an edge in the network,
    a defuault type of INTERACTS will be used.

    Parameters
    ----------
    G : A NetworkX graph
    
    filename : A file location
    
    compressed : Whether to write one edge per line, or to combine edges. (default=False)
    """
    
    with open(filename, 'w') as file:
        if(compressed):
            H = G.copy()
            for node in H.nodes():
                neighbors = H[node]
                types = {}
                for neighbor, attrs in neighbors.iteritems():
                    try:
                        if attrs['type'] in types:
                            types[attrs['type']].append(neighbor)
                        else:
                            types[attrs['type']] = [neighbor]
                    except KeyError:
                        try:
                            types["INTERACTS"].append(neighbor)
                        except KeyError:
                            types["INTERACTS"] = [neighbor]
                for type in types:
                    file.write("%s\t%s\t%s\n" % (node, type, "\t".join(types[type])))
                    for neighbor in types[type]:
                        H.remove_edge(node,neighbor)
        else:
            for edge in G.edges(data=True):
                try:
                    file.write("%s\t%s\t%s\n" % (edge[0],edge[1],edge[2]['type']))
                except KeyError:
                    file.write("%s\tINTERACTS\t%s\n" % (edge[0], edge[1]))
    
###### ADDITIONAL DATA FILE INPUT #######
def readGMTFile(conn, filename):
    """Reads gene sets relating to pathways or otehr annotations from
    a GMT formatted file.

    Note: This method reads a GMT formatted file as downloaded from
    GSEA's annotation portal - which is NOT actual GMT format.

    Paramters
    ---------
    conn : a MySQL connection to the gene database

    filename : a GMT file location

    Returns
    -------
    gmt : a dictionary of pathways IDs to Entrez IDs
    """

    cursor = conn.cursor()
    gmt = {}
    eid_cache = {}
    with open(filename,'r') as file:
        for line in file:
            cols = line.strip().split("\t")

            pathway_name = cols[0]
            pathway_ulr  = cols[1]
            genes        = cols[2:]

            gmt[pathway_name] = set()

            for gene in genes:
                gmt[pathway_name].update(checkEID(cursor,gene))

    return gmt


def readTargetScanFile(conn, filename, species_code=9606):
    """Reads data from a Target Scan output file into memory.

    Parameters
    ----------
    conn : a MySQL connection

    filename : A file location
    
    species_code : Which species to filter for (default=9606).
        Will not filter for species if the value is None.

    Returns
    -------
    TS : A dictionary of lists
        The dictionary key is an miRNA, and its value is a list of targets.
    """
    
    TS = {}
    cursor = conn.cursor()
    p = re.compile('(ENSG\d+)\.?')

    with open(filename) as file:
        next(file) # Skip the header
        for line in file:
            fields = line.strip().split("\t")

            mirna         = fields[0]
            target_id     = fields[1]
            target_symbol = fields[2]
            transcript    = fields[3]
            species       = int(fields[4])
           
            # Skip along if this doesn't meet our filters
            if species_code is not None and species != species_code:
                continue

            eids = getEID(cursor, target_symbol)
            if not eids:
                m = p.match(target_id)
                eids = crossQuery(cursor, 'Ensembl', m.group(1)) 

            if eids:
                if mirna not in TS:
                    TS[mirna] = {}
                TS[mirna].update({eid: True for eid in eids})
    
    # Clean up formatting
    for mirna in TS:
        TS[mirna] = TS[mirna].keys()
    
    return TS
    
###### NETWORK METHODS ########
def getPartners(nodes,CI,depth=1):
    """Returns a list of interacting partners within a defined distance.

    Parameters
    ----------
    nodes : a list of seed nodes to find partners from
    
    CI : the consolidated interactome
    
    depth : how many interactions away we should include (defaut=1)
    
    Returns
    -------
    partner_list : a list containing the input nodes and the partners from the CI up to the defined depth.
    """
    
    neighbors = []
    this_level = []
    next_level = {node: True for node in set(nodes) & set(CI)}
    visited = {}
        
    for i in xrange(depth):
        if not next_level:
            break
        this_level = next_level
        next_level = {}
        while this_level:
            node, trash = this_level.popitem()
            visited[node] = True
            for neighbor in CI[node]:
                if not neighbor in visited:
                    next_level[neighbor] = True
                    
    return visited.keys()

def all_shortest_paths(G, source, target, cutoff=None):
    """Compute all shortest paths in the graph.
    (Modified version of NetworkX's method to include a cutoff value)

    Parameters
    ----------
    G : NetworkX graph

    source : node
        Starting node for path.

    target : node
        Ending node for path.
       
    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    paths: generator of lists
        A generator of all paths between source and target.
    """
    pred = nx.predecessor(G,source,cutoff=cutoff)
    if target not in pred:
        raise nx.NetworkXNoPath()
    stack = [[target,0]]
    top = 0
    while top >= 0:
        node,i = stack[top]
        if node == source:
            yield [p for p,n in reversed(stack[:top+1])]
        if len(pred[node]) > i:
            top += 1
            if top == len(stack):
                stack.append([pred[node][i],0])
            else:
                stack[top] = [pred[node][i],0]
        else:
            stack[top-1][1] += 1
            top -= 1
    
def linkerGraph(startNET,endNET,CI,spLen=None):
    """Merges two graphs and links their nodes along shortest paths.

    Parameters
    ----------
    startNET : The first network to merge
    
    endNET : The second network to merge
    
    CI : The consolidated interactome
    
    spLEN : The maximum shortest path length to include (optional)

    Returns
    -------
    G : A NetworkX Graph
    """
    
    new_nodes = {}

    # iterate through all pairs (a,b), where a is in networkA, b is in networkB
    for start in startNET.nodes(): 
        for end in endNET.nodes(): 
            try:
                paths_gen = all_shortest_path(CI,start,end,cutoff=spLen)
                for path in paths_gen:
                    new_nodes.update({node: True for node in path})
            except nx.NetworkXNoPath:
                continue

    # Get the subgraph of with the shortest path nodes
    SPS = nx.subgraph(CI,new_nodes.keys())
    # Get the nodes from the original graphs
    SPS.add_nodes_from(startNET.nodes())
    SPS.add_nodes_from(endNET.nodes())
    # Get the edges from the original graphs
    SPS.add_edges_from(startNET.edges())
    SPS.add_edges_from(endNET.edges())

    for node in SPS:
        if node in startNET:
            if node in endNET:
                SPS.node[node]['type'] = 'shared'
            else:
                SPS.node[node]['type'] = 'start'
        elif node in endNET:
            SPS.node[node]['type'] = 'end'
        else:
            SPS.node[node]['type'] = 'linker'

    return SPS
    
def pruneGraph(G,CI,goi,pmax=0.05):
    """Prunes a graph of nodes that are not significantly connected to genes of interest.
    Pruning is done using Fisher's Exact Test, comparing the number of interactions a gene has
    to the number of genes of interest it connects to.
    
    This function is a python implementation of NetBox's methods:
    http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008918

    Parameters
    ----------
    G : a graph to prune
    
    CI : the consolidated interactome
    
    goi : a list of genes of interest
    
    pmax : the p-value cutoff for significance (default=0.05)

    Returns
    -------
    P : A NetworkX Graph
    """
    
    import numpy.random as r
    
    goi_set = set(goi)

    N = len(CI) # Population size
    n = len(goi_set & set(CI)) # Successes in population
    
    P = G.copy()

    # Iterate over partner network
    for node in G:

        # Ignore genes of interest
        if node in goi_set:
            continue

        M = len(CI[node]) # Size of sample
        x = len(set(CI[node]) & goi_set) # Successes in sample

        # Compute significance
        p = fishers(N,n,M,x)

        # Prune the node if it is not significant
        if p > pmax:
            P.remove_node(node)

    return P

def pvalGraph(goi,CI,max_p=0.05,max_size=None, verbose=False):
    """Uses the Chan lab's original graph extension method utilizing p-value cutoffs.

    Parameters
    ----------
    goi : a list of genes of interest (GOIs)
    
    CI : the consolidated interactome
    
    max_p : the p-value cutoff for significance (default=0.05)

    max_size : the maximum number of nodes to include in the network (default=2 times the number of GOIs)

    verbose : if True, will print to stdout the progress of the algorithm (default=False)

    Returns
    -------
    G : A NetworkX Graph
    """

    if not max_size:
        max_size = len(goi) * 2

    goi_set = set(goi)
    new_genes = set(goi)
    
    G = CI.subgraph(new_genes)
    comps = list(nx.connected_components(G))

    # Extend while the network is of a reasonable size
    # (Only examines the largest connected component)
    while len(max(nx.connected_components(G), key=len)) < max_size:

        # Consider adding any interacting partner of the current extended network
        neighbors = set(getPartners(new_genes,CI,2)) - new_genes
        # Keep track of the genes with the greatest significance
        best_genes = {'nodes': [],'pval': 2}
        # Our population size is two interactions deep from the current genes
        N = len(getPartners(new_genes,CI,3))

        for neighbor in neighbors:
            # Any included gene must interact with 1 GOI
            # and two genes already in the network
            interactors = set(CI[neighbor])
            if interactors & goi_set == 0 or interactors & new_genes < 2:
                continue
                
            # The gene must be significantly connected to the GOIs
            tp = len(set(CI[neighbor]) & new_genes)
            fp = len(CI[neighbor]) - tp
            fn = len(new_genes) - tp
            tn = N - tp - fp - fn 
            p = pvalue(tp,fp,fn,tn).right_tail
            if p == best_genes['pval']:
                    best_genes['nodes'].append(neighbor)
            elif p < best_genes['pval']:
                best_genes['nodes'] = [neighbor]
                best_genes['pval'] = p
        
        # Stop the process if we fail to add any new nodes
        if len(best_genes['nodes']) == 0 or best_genes['pval'] >= max_p:
            if verbose:
                print "\nTerminating: All significant nodes included"
            break

        new_genes.update(best_genes['nodes'])
        G = CI.subgraph(new_genes)
        if verbose:
            ccs = list(nx.connected_components(G))
            LCC = max(ccs, key=len)
            sys.stdout.write("Number of Components: %d\tLCC size: %d\t%d/%d seeds included\r" % (len(ccs),len(LCC),len(set(LCC)& goi_set),len(goi)))
            sys.stdout.flush()
       

    if verbose:
        if len(max(nx.connected_components(G), key=len)) >= max_size:
            print "\nTerminating: Maximum network size reached"
        else:
            print


    # Return only the largest connected component
    return CI.subgraph(max(nx.connected_components(G), key=len))


def spGraph(goi, CI, iters=1000, prune=True, verbose=False):
    """Creates a pruned, expanded network from the spNodes function.
    This method takes the nodes from spNodes, given the goi list and
    the CI graph, and removes extraneous nodes and edges. Using Monte
    Carlo estimation, p-values are assigned to all edges in the graph.
    Edges that are not statistically significant are then pruned to
    reduce potential noise and to make clustering easier.

    Paramters
    ---------
    goi : a list of genes of interest

    CI : the consolidated interactome

    iters : the number of iterations to use for Monte Carlo estimation.
            A large number of iterations increases accuracy of p-values,
            but also takes longer to run. (default=1000)

    prune : by default, this method will use the Monte Carlo method to
            prune the graph. However, if you just want a subgraph of
            the nodes from spNodes, change this parameter to False.
            (default=True)

    verbose : if True, will print progress of Monte Carlo estimation
              to STDERR. (default=False)

    Returns
    -------
    G : an expanded NetworkX Graph
    """
 
    # The maximum frequency we can have
    #max_freq = iters * max_p

    # Create the spGraph for the GOIs
    G = CI.subgraph(spNodes(goi,CI))
    if not prune:
        return G

    # Use Monte Carlo estimation to prune edges
    genes = CI.nodes()
    n = len(goi)
    for i in xrange(iters):
        if verbose:
            sys.stderr.write("%d/%d\r" % (i+1, iters))
            sys.stderr.flush()
        # Create an spGraph from random genes
        random.shuffle(genes)
        O = spNodes(genes[:n],CI)

        #for edge in O.edges():
        for edge in G.subgraph(O).edges():
            # Count the frequency by which each edge
            # in G appears at random
            node1, node2 = edge
            if 'pval' in G[node1][node2]:
                G[node1][node2]['pval'] += 1
                # If an edge appears too pvaluently, prune it
                #if G[node1][node2]['pval'] >= max_pval:
                #    G.remove_edge(*edge)
            else:
                G[node1][node2]['pval'] = 1
    if verbose:
        sys.stderr.write("\n")

    # Calculate p-values
    edges = G.edges(data=True)
    for edge in edges:
        if 'pval' in edge[2]:
            edge[2]['pval'] = (1 + edge[2]['pval']) / float(iters + 1)
        else:
            edge[2]['pval'] = 1 / float(iters + 1)

    # Start pruning
    edges.sort(key=lambda x: x[2]['pval'], reverse=True)
    for edge in edges:
        # Stop if we've pruned all insignificant edges
        if edge[2]['pval'] < 0.05:
            break

        # Check to see if pruning this edge makes a GOI orphan
        if edge[0] in goi and len(G[edge[0]]) == 1:
            break
        if edge[1] in goi and len(G[edge[1]]) == 1:
            break

        # Otherwise, remove the insignificant edge
        G.remove_edge(*edge[:2])

    return max(nx.connected_component_subgraphs(G), key=len)

def spNodes(goi, CI):
    """Generates a list of nodes expanded from genes of interest (GOI). The
    process starts by selecting the largest group of GOI that interact with
    each other. From there, each non-connected GOI attempts to connect using
    at most 1 linker genes (i.e. GOI <-> Linker <-> GOI). This process is
    repeated until either all the GOI are included or there are no more GOI
    that are only 1 linker gene away.

    Parameters
    ----------
    goi : a list of genes of interest

    CI :  the consolidated interactome to take interactions from

    Returns
    -------
    nodes : an expanded list of nodes related to the genes of interest
    """

    # Nodes to return
    nodes = set(max(nx.connected_components(CI.subgraph(goi)), key=len))

    # Keep track of which targets need to be visited
    this_level = None
    next_level = nodes.copy()

    # Keep track of orphans
    orphans = set(goi) - next_level

    # Iterate until we run out of orphans, or we run out of connections
    while next_level and next_level:
        this_level = next_level
        next_level = set()

        # Try to link each orphan to connected nodes
        for orphan in orphans:

            # To decrease lookup time, store the partners now
            partners = set(CI[orphan])

            # Try to connect orphan to any/all nodes in the LCC
            # (So long as we haven't already tried)
            for target in this_level:
                linkers = partners & set(CI[target])
                # If we made a successful link
                if linkers:
                    # Add the linker genes and orphan to the next leve
                    next_level.update(linkers)
                    next_level.add(orphan)
                    # And flag the linker to be removed

        orphans -= next_level
        nodes.update(next_level)

    return nodes


def infomapCluster(G, weight_feature=None):
    """Computes clusters using the infomap algorithm

    Parameters
    ----------
    G : A NetworkX Graph to cluster
        This object will also be modified. Nodes will gain
        a new attribute called 'cluster' indicating which
        cluster it belongs to.

    weight_feature : If given, infomap will add edge weights
                     from the given feature found in G.
                     (default=None)

    Returns
    -------
    clustDict : A dictionary of clusters to node lists
    """

    idToNode = G.nodes()
    nodeToId = {idToNode[i] : i for i in xrange(len(idToNode))}

    fname = uuid4()
    with open('/tmp/%s.llf' % (fname),'w') as file:
        for edge in G.edges(data=True):
            if weight_feature is None:
                file.write("%d %d 1\n" % (nodeToId[edge[0]], nodeToId[edge[1]]))
            else:
                file.write("%d %d %d\n" % (nodeToId[edge[0]], nodeToId[edge[1]], edge[2][weight_feature]))

    proc = subprocess.Popen([
            "/home/HandenA/extend_network/Infomap/Infomap.exe",
            "--input-format", "link-list",
            "--zero-based-numbering",
            "--clu",
            "--undirected",
            "--silent",
            "/tmp/%s.llf" % (fname),
            "/tmp/"]).wait()

    clustDict = {}
    with open('/tmp/%s.clu' % (fname),'r') as file:
        for line in file:
            if line[0] == "#":
                continue
            id, cluster, flow  = line.strip().split(" ")
            id = int(id)
            cluster = int(cluster)
            G.node[idToNode[id]]['cluster'] = cluster
            if cluster not in clustDict:
                clustDict[cluster] = []
            clustDict[cluster].append(idToNode[id])

    os.remove('/tmp/%s.llf' % (fname))
    os.remove('/tmp/%s.clu' % (fname))

    return clustDict

def louvainCluster(G):
    """Computes clusters using the Louvain algorithm

    Parameters
    ----------
    G : A NetworkX Graph to cluster
        This object will also be modified. Nodes will gain
        a new attribute called 'cluster' indicating which
        cluster it belongs to.

    Returns
    -------
    clustDict : A dictionary of clusters to node lists
    """

    try:
        clust = community.best_partition(G)    # attempt louvain method
    except:
        clust = {}                # if clustering fails, assign all nodes to the same cluster
        for x in G:
            clust[x] = 0

    clustDict = {}
    for x in clust:
        G.node[x]['cluster'] = clust[x] # tag nodes by clusterID
        if clust[x] in clustDict: # rework dictionary
            clustDict[clust[x]].append(x)
        else:
            clustDict[clust[x]] = [x]

    return clustDict

##### STATISTIC METHODS ######
def fishers(N,n,M,x):
    """Estimates the p-value of Fisher's Exact Test

    Parameters
    ----------
    N : Population size
    
    n : Number of successes in the population
    
    M : Sample size
    
    x : Number of successes in the sample
    
    Returns
    -------
    p : the right-tailed p-value
    """
    tp = x
    fp = M - x
    fn = n - x
    tn = N - tp - fp - fn
    return pvalue(tp,fp,fn,tn).right_tail
    
def spanningScores(TS, NET, CI, clusters):
    """Computes spanning scores for all miRNA.

    Parameters
    ----------
    TS : the TargetScan dictionary

    NET : the network to compute scores over

    CI : the consolidated interactome

    clusters : a list of clusters for NET

    Returns
    -------
    """

    stats = {}
    N = len(CI)  # population size
    M = len(NET) # sample size

    data = {}
    for miR, targets in TS.iteritems():
        targpool = set(targets) & set(CI) # get target pool

        n = len(targpool)                # successes in population
        x = len(targpool & set(NET))    # successes in sample

        pval = fishers(N,n,M,x)
        negp = -log(pval)
        
        clust_count = sum([1 for cluster in clusters if len(cluster) >3 and set(targets) & set(cluster)])

        data[miR] = {
            'hits'        : x,
            'negp'        : negp,
            'clust_count' : clust_count,
            'pval'        : pval
        }

    max_o = float(max([data[miR]['hits']        for miR in data])) # best achieved overlap
    max_c = float(max([data[miR]['clust_count'] for miR in data])) # best achieved cluster count
    max_p = float(max([data[miR]['negp']        for miR in data])) # best achieved -log(pval)

    if max_o == 0:
        raise Exception("There are no miRNA targets in the network")
    if max_c == 0:
        raise Exception("There are no miRNA targets in any of the clusters")

    scores = {}
    for miR, stats in data.iteritems():

        first = stats['negp'] / max_p    # -log(pval) / -log(best_pval) (max = 1)
        second = stats['clust_count'] / max_c # cluster_count / best_cluster_count (max = 1)
        third = stats['hits'] / max_o # hit_count / best_hit_count (max = 1)

        scores[miR] = first + second + third    # sum score components (max score = 3)

    return scores
