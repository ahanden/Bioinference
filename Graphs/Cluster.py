from subprocess import Popen
from uuid import uuid4
import community
import os

def infomap(graph, weight_feature=None):
    """Computes clusters using the infomap algorithm

    Parameters
    ----------
    graph : A NetworkX Graph to cluster. This object will also be modified.
            Nodes will gain a new attribute called 'cluster' indicating which
           cluster it belongs to.

    weight_feature : If given, infomap will add edge weights
                     from the given feature found in G.
                     (default=None)

    Returns
    -------
    A dictionary of clusters to node lists
    """

    idToNode = list(graph)
    nodeToId = {idToNode[i] : i for i in xrange(len(idToNode))}

    fname = uuid4()
    with open('%s.llf' % (fname),'w') as file:
        for edge in graph.edges(data=True):
            if weight_feature is None:
                file.write("%d %d 1\n" % (nodeToId[edge[0]], nodeToId[edge[1]]))
            else:
                file.write("%d %d %d\n" % (nodeToId[edge[0]], nodeToId[edge[1]], edge[2][weight_feature]))

    proc = Popen([
            "C:\\infomap\\Infomap.exe",
            "--input-format", "link-list",
            "--zero-based-numbering",
            "--clu",
            "--undirected",
            "--silent",
            "%s.llf" % (fname),
            "./"]).wait()

    clustDict = {}
    with open('%s.clu' % (fname),'r') as file:
        for line in file:
            if line[0] == "#":
                continue
            id, cluster, flow  = line.strip().split(" ")
            id = int(id)
            cluster = int(cluster)
            graph.node[idToNode[id]]['cluster'] = cluster
            if cluster not in clustDict:
                clustDict[cluster] = []
            clustDict[cluster].append(idToNode[id])

    os.remove('%s.llf' % (fname))
    os.remove('%s.clu' % (fname))

    return clustDict

def louvain(graph):
    """Computes clusters using the Louvain algorithm

    Parameters
    ----------
    graph : A NetworkX Graph to cluster. This object will also be modified.
            Nodes will gain a new attribute called 'cluster' indicating which
            cluster it belongs to.

    Returns
    -------
    A dictionary of clusters to node lists
    """

    try:
        clust = community.best_partition(graph)    # attempt louvain method
    except:
        clust = {}                # if clustering fails, assign all nodes to the same cluster
        for x in graph:
            clust[x] = 0

    clustDict = {}
    for x in clust:
        graph.node[x]['cluster'] = clust[x] # tag nodes by clusterID
        if clust[x] in clustDict: # rework dictionary
            clustDict[clust[x]].append(x)
        else:
            clustDict[clust[x]] = [x]

    return clustDict
