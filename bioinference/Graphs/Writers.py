#!/bin/python

import networkx as nx

def write_tsv(stream, G, delimiter="\t"):
    """Writes a network's edges to a stream

    Parameters
    ----------
    stream : Any object with a write() method

    G : A NetworkX Graph

    delimiter : The file delimiter (default=tab)

    """

    for edge in G.edges():
        stream.write("%s%s%s\n" % (edge[0], delimiter, edge[1]))

def write_sif(stream, G):
    """Writes a Graph to a stream in SIF format

    Parameters
    ----------
    stream : Any object with a write() method

    G : A NetworkX graph
    """
    
    for edge in G.edges(data=True):
        typ=edge[2].get("type", "INTERACTS")
        stream.write("%s\t%s\t%s\n" % (edge[0], typ, edge[1]))
