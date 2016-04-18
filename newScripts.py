import networkx as nx
from fisher import pvalue
import MySQLdb

####### NETWORK LOADING FUNCTIONS #########
def readDB(conn):
    """Reads interactions into a NetworkX graph from a MySQL database

    Parameters
    ----------
    conn : A MySQLdb connection

    Returns
    -------
    G : A NetworkX Graph
    """
    G = nx.Graph()
    cursor = conn.cursor()
    cursor.execute("SELECT entrez_id1, entrez_id2 FROM interactions")
    for row in cursor.fetchall():
        G.add_edge(row[0],row[1])
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
        node1    node1    EntrezID1    EntrezID2
    The file should NOT have a header line.

    Parameters
    ----------
    filename : A file location

    Returns
    -------
    G : A NetworkX Graph
    """
    
    G = nx.Graph()
    with open(filename,'r') as file:
        for line in file:
            (idA, idB, nodeA, nodeB) = line.strip().split("\t")
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

def extendGraph(goi,CI,max_p=0.05,max_size=None):
    """Uses the Chan lab's original graph extension method utilizing p-value cutoffs.

    Parameters
    ----------
    goi : a list of genes of interest (GOIs)
    
    CI : the consolidated interactome
    
    max_p : the p-value cutoff for significance (default=0.05)

    max_size : the maximum number of nodes to include in the network (default=2 times the number of GOIs)

    Returns
    -------
    G : A NetworkX Graph

    TODO
    ----
    Decide whether to add nodes one at a time or as a bunch
    Increase efficiency by skipping obviously insignificant nodes
    """

    if not max_size:
        max_size = len(goi) * 2

    goi_set = set(goi)
    new_genes = set(goi)
    
    G = CI.subgraph(new_genes)
    comps = list(nx.connected_components(G))

    # Extend while the network is of a reasonable size
    # (Only examines the largest connected component)
    print "Starting with a LCC of size %d" % (len(max(nx.connected_components(G), key=len)))
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
            print "Maxed out significant genes"
            break

        new_genes.update(best_genes['nodes'])
        G = CI.subgraph(new_genes)
        print "LCC has grown to %d nodes" % (len(max(nx.connected_components(G), key=len)))
        
    # Return only the largest connected component
    return CI.subgraph(max(nx.connected_components(G), key=len))

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
