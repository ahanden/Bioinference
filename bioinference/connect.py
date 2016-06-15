import networkx as nx
import random

def random_walk(G, genes, pr=0.333, iters=1000, always_walk=True):
    """Ranks the nodes in G based on their probability of
    being visited on random walks.

    Parameters
    ----------
    G - The Networkx Graph to traverse

    genes - a list of nodes to start from 

    pr - the probability of restart (default=0.333)

    iters - the number of iterations to do for each node in genes (default=1000)

    always_walk - whether the walk should traverse an edge on ever time step or
                  be able to stay in one place at random (default=True)

    Returns
    -------
    A dictionary of nodes to probabilities of being visited
    """

    ranked_nodes = {n: 0 for n in G}

    for gene in set(genes):
        del ranked_nodes[gene]
        node = gene
        for i in xrange(iters):
            if always_walk:
                node = random.choice(G[node].keys()) if random.random() > pr else gene
            else
                node = random.choice(G[node].keys() + [node]) if random.random() > pr else gene

            if node in ranked_nodes:
                ranked_nodes[node] += 1

    return {n: v/float(len(set(genes)) * iters) for n, v in ranked_nodes.iteritems()}
