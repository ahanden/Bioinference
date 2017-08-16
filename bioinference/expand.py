from fisher import pvalue
from random import shuffle
import networkx as nx
import sys
from itertools import chain

def prunedGraph(goi, CI, pmax=0.05, G=None):
    """Creates am expanded graph of nodes that are significantly connected to genes of
    interest. 
    
    The algorithm creates an initial graph of all genes of interest and their immediate 
    interactors. Then, the immediate interactors are pruned based on Fisher's Exact Test, 
    comparing the number of interactions a gene has to the number of genes of interest it 
    connects to. Genes with a p-value less than pmax are removed from the graph.
    
    This function is a python implementation of NetBox's methods.

    Cerami E, Demir E, Schultz N, Taylor BS, Sander C. Automated network analysis 
    identifies core pathways in glioblastoma. PLoS One. 2010 Feb 12;5(2):e8918. doi: 
    10.1371/journal.pone.0008918. PubMed PMID: 20169195; PubMed Central PMCID:
    PMC2820542.
    (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008918)

    Parameters
    ----------
    goi : a list of genes of interest
    
    CI : the consolidated interactome
    
    pmax : the p-value cutoff for significance (default=0.05)

    G : you may supply a pre-generated graph to prune nodes from instead of
        using this function's way of generating the initial graph (optional)

    Returns
    -------
    P : A NetworkX Graph
    """
    
    if G is None:
        G = CI.subgraph(getPartners(goi,2))

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

def cottrillGraph(goi, CI, max_p=0.05, max_size=None, verbose=False):
    """Creates an expanded graph by iteratively adding genes based on the
    results of a Fisher's Exact Test.

    This method starts with a graph of all the given genes of interest. 
    All immediate interactors of the genes of interest are then ranked 
    based on the results of Fisher's Exact Test - comparing the number 
    of genes of interest a gene interacts with to the total number of 
    interactions that gene has. The most significant gene(s) are added
    to the graph. 
    
    The process is then repeated, treating the added gene(s) as genes of 
    interest as well. The expansion terminates when either there are no
    more significant genes left (p > max_p) or the largest connected 
    component of the network has grown too large (|G| >= max_size).

    Only the largest connected component is returned.

    This is a re-write of the methods used by Alex Cottrill.

    Parikh VN, Jin RC, Rabello S, Gulbahce N, White K, Hale A, Cottrill KA, Shaik 
    RS, Waxman AB, Zhang YY, Maron BA, Hartner JC, Fujiwara Y, Orkin SH, Haley KJ,
    Barabasi AL, Loscalzo J, Chan SY. MicroRNA-21 integrates pathogenic signaling to 
    control pulmonary hypertension: results of a network bioinformatics approach.
    Circulation. 2012 Mar 27;125(12):1520-32. doi: 10.1161/CIRCULATIONAHA.111.060269.
    Epub 2012 Feb 27. PubMed PMID: 22371328; PubMed Central PMCID: PMC3353408.
    (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3353408/)

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


    TODO
    ----
    Think of a good name for this method
    """

    if max_size is None:
        max_size = len(goi) * 2

    goi_set = set(goi)
    new_genes = goi_set & set(CI)
    
    G = CI.subgraph(new_genes)

    # Extend while the network is of a reasonable size
    # (Only examines the largest connected component)
    while len(max(nx.connected_components(G), key=len)) < max_size:

        # Consider adding any interacting partner of the current extended network
        neighbors = set(chain(*list(CI[g].keys() for g in new_genes))) - new_genes

        # Keep track of the genes with the greatest significance
        best_genes = {'nodes': [],'pval': 1}

        # Our population size is two interactions deep from the current genes
        #N = len(getPartners(new_genes,CI,3))
        N = len(neighbors | new_genes | set(chain(*list(CI[g].keys() for g in neighbors))))

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
            sys.stdout.write("Number of Components: %d, LCC size: %d, %d/%d seeds included     \r" % (len(ccs),len(LCC),len(set(LCC)& goi_set),len(goi)))
            sys.stdout.flush()
       

    if verbose:
        if len(max(nx.connected_components(G), key=len)) >= max_size:
            print "\nTerminating: Maximum network size reached"
        else:
            print


    # Return only the largest connected component
    return CI.subgraph(max(nx.connected_components(G), key=len))


def linkerGraph(goi, CI):
    """Creates an expanded graph by iteratively connecting genes of
    interest into a largest connected component (LCC).

    An initial graph is created by finding the largest connected
    component of genes of interest only. Each unincluded GOI is
    included by finding any intermediate interactor to connect
    it to the LCC. That is, if there is path such that
    Unincluded GOI <-> Other Gene <-> Included GOI, the path is
    added to the network. This process is repeated until either
    all GOIs are included, or no more GOIs can be included.
    
    Paramters
    ---------
    goi : a list of genes of interest

    CI : the consolidated interactome

    Returns
    -------
    G : an expanded NetworkX Graph
    """

    included = set(max(nx.connected_components(CI.subgraph(goi)), key=len))
    orphaned = set(goi) - included

    # Keep track of which targets need to be visited
    this_level = None
    next_level = included.copy()

    # Iterate until we run out of orphans, or we run out of connections
    while next_level and orphaned:

        this_level = next_level
        next_level = set()

        # Try to link each orphan to connected nodes
        for orphan in orphaned:

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

        orphaned -= next_level
        included.update(next_level)

    return CI.subgraph(included)


def monteCarloPrune(goi, CI, method, iters=1000, max_p=0.05, *args):
    """Creates a pruned graph by using Monte Carlo estimation.

    This method first generates an expanded graph with the user-
    provided method. It will then generate a defined number of random
    networks using the same method and random genes of interest equal
    in number to the original genes of interest.

    Each edge in the original expanded graph is then given a p value
    based on how frequently it appears in the random graphs. The
    edges are sorted in descending order by p value and are pruned.
    Edges continue to be pruned until either all remaining edges are
    lower than max_p or removing another edge would orphan a GOI. 

    Paramters
    ---------
    goi : a list of genes of interest

    CI : the consolidated interactome

    method : the graph expansion method to use for generating graphs

    iters : the number of iterations to use for Monte Carlo estimation.
            A large number of iterations increases accuracy of p-values,
            but also takes longer to run. (default=1000)

    max_p : the threshold at which to stop pruning (default=0.05)

    *args : arguments to be passed to the graph expansion method

    Returns
    -------
    G : an expanded, pruned NetworkX Graph
    """
 
    goi = set(goi)

    # Create the spGraph for the GOIs
    G = method(goi,CI,*args)
    for edge in G.edges(data=True):
        edge[2]['pval'] = 0

    # Use Monte Carlo estimation to prune edges
    genes = CI.nodes()
    n = len(goi)
    for i in xrange(iters):
        # Create a graph from random genes
        shuffle(genes)
        O = method(genes[:n],CI,*args)

        # Count the frequency by which each edge
        # in G appears at random
        for edge in O.edges():
            if G.has_edge(*edge):
                node1, node2 = edge
                if 'pval' in G[node1][node2]:
                    G[node1][node2]['pval'] += 1
                else:
                    G[node1][node2]['pval'] = 1

    # Start pruning
    edges = G.edges(data=True)
    edges.sort(key=lambda x: x[2]['pval'], reverse=True)
    max_freq = (1 + iters) * max_p - 1
    for edge in edges:
        # Stop if we've pruned all insignificant edges
        if edge[2]['pval'] < max_freq:
            break

        # Check to see if pruning this edge makes a GOI orphan
        T = G.copy()
        T.remove_edge(*edge[:2])
        T = next(nx.connected_component_subgraphs(T))

        if len(set(T) & goi) < len(set(G) & goi):
            break

        # Otherwise, remove the insignificant edge
        G = T

    # Calculate p-values
    edges = G.edges(data=True)
    for edge in edges:
        if 'pval' in edge[2]:
            edge[2]['pval'] = (1 + edge[2]['pval']) / float(iters + 1)
        else:
            edge[2]['pval'] = 1 / float(iters + 1)

    return max(nx.connected_component_subgraphs(G), key=len)


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

