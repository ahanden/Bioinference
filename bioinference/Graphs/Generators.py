#!/usr/bin/python

from fisher import pvalue
from random import shuffle
import networkx as nx

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
    The right-tailed p-value
    """
    
    tp = x
    fp = M - x
    fn = n - x
    tn = N - tp - fp - fn
    return pvalue(tp,fp,fn,tn).right_tail

def pruned_graph(seeds, interactome, p_max=0.05, graph=None):
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
    seeds : An iterable of Genes to construct the graph from
   
    interactome : A NetworkX Graph of the interactome
    
    p_max : the p-value cutoff for significance (default=0.05)

    graph : you may supply a pre-generated graph to prune nodes from instead of
            using this function's way of generating the initial graph (optional)

    Returns
    -------
    A NetworkX Graph
    """

    seeds = set(seeds)
   
    if graph is None:
        nodes = seeds
        for seed in seeds:
            nodes.update(interactome[seed].keys())
        graph = interactome.subgraph(nodes)

    N = len(interactome) # Population size
    n = len(seeds & set(interactome)) # Successes in population

    # Iterate over partner network
    for node in nodes - seeds:

        M = len(interactome[node]) # Size of sample
        x = len(set(interactome[node]) & seeds) # Successes in sample

        # Compute significance
        p = fishers(N,n,M,x)

        # Prune the node if it is not significant
        if p >= pmax:
            graph.remove_node(node)

    return graph

def cottrill_graph(seeds, interactome, max_p=0.05, max_size=None):
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
    seeds : An iterable of genes to generate the graph with
    
    interactome : The consolidated interactome
    
    max_p : The p-value cutoff for significance (default=0.05)

    max_size : The maximum number of nodes to include in the network (default=2 times the number of GOIs)

    Returns
    -------
    A NetworkX Graph
    """

    seeds = set(seeds)

    if max_size is None:
        max_size = len(seeds) * 2

    new_genes = seeds & set(CI)
    
    graph = interactome.subgraph(new_genes)

    # Extend while the network is of a reasonable size
    # (Only examines the largest connected component)
    while len(max(nx.connected_components(graph), key=len)) < max_size:

        # Find all immedaite interactors of the network in the CI taht are
        # not already included
        neighbors = set()
        for node in graph:
            neighbors.update(interactome[node])
        neighbors -= new_genes

        # Keep track of the genes with the greatest significance
        best_genes = {'nodes': [],'pval': 1}

        # Our population size is two interactions deep from the current genes
        #N = len(getPartners(new_genes,CI,3))
        N = set()
        for node in neighbors:
            N.update(interactome[node])
        N = len( N | neighbors | new_genes )

        for neighbor in neighbors:
            # Any included gene must interact with 1 seed and 2 genes 
            # already in the network
            interactors = set(CI[neighbor])
            if interactors & seeds == 0 or interactors & new_genes < 2:
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
            break

        new_genes.update(best_genes['nodes'])
        graph = interactome.subgraph(new_genes)

    if verbose:
        if len(max(nx.connected_components(G), key=len)) >= max_size:
            print "\nTerminating: Maximum network size reached"
        else:
            print


    # Return only the largest connected component
    return CI.subgraph(max(nx.connected_components(G), key=len))


def linker_graph(seeds, interactome):
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
    seeds : An iterable of genes to generate the network from

    interactome : The consolidate interactome

    Returns
    -------
    An expanded NetworkX Graph
    """

    seeds = set(seeds)

    nodes = set(max(
        nx.connected_components(
            interactome.subgraph(seeds)),
        key=len))

    orphaned = seeds - nodes

    # Keep track of which targets need to be visited
    this_level = None
    next_level = included.copy()

    while orphaned:
        # Find any common interactors for the orphans and the current graph
        linkers = set(interactome[orphaned]) & set(interactome[nodes])

        # Give up if we can't find any new connections
        if not linkers:
            break

        # Update the graph to include the old graph, new linkers and any
        # newly included orphans
        nodes = set(max(
            nx.connected_components(
                interactome.subgraph(nodes | linkers | orphans)),
            key=len))

        orphaned -= nodes

    return interactome.subgraph(nodes)


def monte_carlo_pruned_graph(seeds, interactome, method, iters=1000, max_p=0.05, *args):
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

    NOTE: This method requires that whatever method you use to create the
          initial graph, it must be pruned to a largest connected component.

    Paramters
    ---------
    seeds : An iterable of genes to generate the network from

    interactome : The consolidate interactome

    method : The graph expansion method to use for generating graphs

    iters : The number of iterations to use for Monte Carlo estimation.
            A large number of iterations increases accuracy of p-values,
            but also takes longer to run. (default=1000)

    max_p : The threshold at which to stop pruning (default=0.05)

    *args : Arguments to be passed to the graph expansion method

    Returns
    -------
    An expanded, pruned NetworkX Graph
    """

    seeds = set(seeds)

    # Create the initial Graph, all p-values assumed to be minimum
    graph = max(
        nx.connected_component_subgraphs(
            method(seeds, interactome, *args)),
        key=len)
    for edge in graph.edges(data=True):
        edge[2]['pval'] = 1.0 / (iters + 1)

    # Use Monte Carlo estimation to prune edges
    genes = interactome.nodes()
    n = len(seeds)
    for i in xrange(iters):
        # Create a graph from random genes
        shuffle(genes)
        trial = max(
            nx.connected_component_subgraphs(
                method(genes[:n], interactome, *args)),
            key=len)

        # Count the frequency by which each edge
        # in the graph appears at random
        for edge in trial.edges():
            if graph.has_edge(*edge):
                node1, node2 = edge
                if 'pval' in graph[node1][node2]:
                    graph[node1][node2]['pval'] += 1.0 / (iters + 1)

    # Start pruning
    for edge in sorted(
            graph.edges(data=True),
            key=lambda x: x[2]['pval'],
            reverse=True):

        # Skip if we already pruned this edge
        if not graph.has_edge(*edge[:2]):
            continue

        # Stop if we've pruned all insignificant edges
        if edge[2]['pval'] < max_p:
            break


        # Check to see if pruning this edge makes a GOI orphan
        trial = graph.copy()
        trail.remove_edge(*edge[:2])
        trail = max(nx.connected_component_subgraphs(T), key=len)

        if len(set(trial) & seeds) < len(se(graph) & seeds):
            break

        # Otherwise, remove the insignificant edge
        graph = trial

    return graph
