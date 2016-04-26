from fisher import pvalue
from subprocess import call
from uuid import uuid4
import MySQLdb
import community
import networkx as nx
import os
import sys

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
    #cursor.execute("SELECT entrez_id1, entrez_id2, count(pubmed_id) FROM interactions LEFT JOIN publications ON interactions.int_id = publications.int_id GROUP BY interactions.int_id")
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
        
        line = file[0]    
        if "\t" in line:
            delimiter = "\t"
            
        for line in file:
            fields = line.strip().split(delimiter)
            
            if len(fields) < 3:
                raise IndexError('The SIF file is improperly formatted. Visit http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats')
            
            eid1 = fields.pop(0)
            type = fields.pop(0)
            
            for eid2 in fields:
                if eid1 != eid2:
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
                    file.write("%s\tINTERACTSt%s\n" % (edge[0], edge[1]))
    
###### ADDITIONAL DATA FILE INPUT #######
def readTargetScanFile(filename, species_code="9606"):
    """Reads data from a Target Scan output file into memory.

    Parameters
    ----------
    filename : A file location
    
    species_code : Which species to filter for (default=9606).
        Will not filter for species if the value is None.

    Returns
    -------
    TS : A dictionary of lists
        The dictionary key is an miRNA, and its value is a list of targets.
    """
    
    TS = {}
    species_code = str(species_code)

    with open(filename) as file:
        for line in file:
            (mirna, garbage, target, garbage, species) = line.strip().split("\t")
            
            if species_code == "None" or species == species_code:
                if mirna in TS:
                    TS[mirna][target] = True
                else:
                    TS[mirna] = {target: True}
    
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
        
    for i in range(depth):
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

def extendGraph(goi,CI,max_p=0.05,max_size=None, verbose=False):
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
       

    if verbose and len(max(nx.connected_components(G), key=len)) >= max_size:
        print "\nTerminating: Maximum network size reached"
    elif verbose:
        print


    # Return only the largest connected component
    return CI.subgraph(max(nx.connected_components(G), key=len))

def spGraph(goi, CI, max_dist=None, min_pubs=0, filter=False):
    """Creates a network from genes of interest (GOIs) based on shortest paths.
    The process begins by creating a subnetwork of interactions just among GOIS.
    The largest connected component (LCC) of this network is used as the base.
    Then, for as many iterations as defined, orphaned GOIs will be included that
    are exactly one interaction away from the LCC.

    Parameters
    ----------
    goi : a list of genes of interest (GOIs)
    
    CI : the consolidated interactome

    max_dist : If an int, the number of iterations to run the expansion process.
               If None, the algorithm will continue to expand until it cannot
               include any more GOIs. (default=None, min=1)

    min_pubs : The minimum number of publications an edge must have to be added to the
               LCC. Note that this is only for expanding the network, and does not
               affect the initial LCC made from all the GOIs. (default=0)

    filter : If True, only the edges with the most publications will be used for
             expansion.
             If False, all edges connecting orphan nodes will be used. (default=False)
                

    Returns
    -------
    G : A NetworkX Graph of the expanded network
    """

    # Check that max_dist is acceptable
    if max_dist is not None and max_dist < 1:
        raise Exception('max_dist must be a positive integer or None.')

    # Graph structure to eventually be returned
    G = nx.Graph()

    # Genes to seed with for each iteration
    seeds = set(goi)

    # Keep track of how many seeds we had in the last iteration
    seed_tracker = 0

    # Keep track of how many iterations we've had
    iter_count = 0

    # Keep the loop while the following 3 conditions are *all* met
    #   - The maximum distance has not been reached (or no maximum distance was provided)
    #   - The expansion algorithm is still adding new genes
    #   - There are still orphan nodes left
    while (maxDist is None or iter_count < maxDist) and seed_tracker < len(seeds) and len(set(goi) - set(G)) > 0:
        # Keep track of our place, in case a termination value is given
        iter_count += 1

        # Keep track of the number of seeds
        seed_tracker = len(seeds)

        # Start by creating a subgraph made of *just* seeds
        # and only use the LCC
        G = max(nx.connected_component_subgraphs(CI.subgraph(goi), key=len))

        # Identify the orphaned GOIs
        orphans = set(goi) - set(G)
        # Also keep track of the included GOIs
        included = set(G)

        # Try to include each orphan...
        for source in orphans:
            # ... by finding paths to any gene already in the network
            for target in included:
                # If we are only including the most confident expandable edges
                if filter:
                    best_linkers = []
                    most_pubs = min_pubs
                    for linker in set(CI[source]) & set(CI[target]):
                        pubs = CI[source][linker]['publications']
                        if pubs > most_pubs:
                            best_linkers = [linker]
                        elif pubs == most_pubs:
                            best_linkers.append(linker)
                    seeds.update(best_linkers) 
                # If we are including all expandable edges
                else:
                    seeds.update(set(CI[source]) & set(CI[target]))

    return G

def infomapCluster(G):
    """Computes clusters using the infomap algorithm

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

    idToNode = G.nodes()
    nodeToId = {idToNode[i] : i for i in range(len(idToNode))}

    fname = uuid4()
    with open('/tmp/%s.llf' % (fname),'w') as file:
        for edge in G.edges():
            file.write("%d %d 1\n" % (nodeToId[edge[0]], nodeToId[edge[1]]))

    call("Infomap/Infomap --input-format link-list --zero-based-numbering --clu --undirected --silent /tmp/%s.llf /tmp/" % (fname),shell=True)

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
    
def pValMIR(miR,NET,TS,CI):
    """Finds the TargetScan targets for a given miRNA found in a given
    network and the consolidated interactome.

    Parameters
    ----------
    miR : miRNA family
    
    NET : the network
    
    TS : TargetScan dictionary
    
    CI : The consolidated interactome
    
    Returns
    -------
    p : p-value computed from a hypergeometric distribution
    
    x : the size of the overlap
    
    TODO
    ----
    Switch to more efficient hypergeometric library.
    """

    targpool = set(TS[miR]) & set(CI)

    N = len(CI)                        # population size
    n = len(targpool)                # successes in population
    M = len(NET)                    # sample size
    x = len(targpool & set(NET))    # successes in sample

    p = fishers(N,n,M,x)

    return (p,x)
    
#def
