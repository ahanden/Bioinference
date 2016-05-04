#!/bin/python

import MySQLdb
import networkx as nx
import genes

"""
Methods for reading and writing graph structures
"""


def mapIDs(G, conn, in_type="Entrez", out_type="Entrez"):
    """Creates a version of the given graph with mapped gene identifiers.

    Parameters
    ----------
    G : a NetworkX graph to translate

    conn : a MySQLdb connection to the genes database

    in_type : the type of ID current used for the graph (default="Entrez")

    out_type : the type of ID to use for the returned graph

    Returns
    -------
    H : a copy of G with out_type identifiers for nodes
    """

    # Map all of the IDs
    id_map = {}
    for node in G:
        eids = []

        if in_type == "Entrez":
            eids = genes.checkEID(node, conn)
        elif in_type == "Symbol":
            eids = genes.getEID(node, conn)
        else:
            eids = genes.crossQuery(node, in_type, conn)

        if out_type == "Entrez":
            id_map[node] = eids
        elif out_type == "Symbol":
            symbols = set()
            for eid in eids:
                symbols.update(genes.getSymbol(eid, conn))
            id_map[node] = symbols
        else:
            fids = set()
            for eid in eids:
                fids.update(genes.crossQuery(node, out_type, conn))
            id_map[node] = fids

    # Re-construct the graph
    H = nx.Graph()
    for edge in G.edges(data=True):
        (node1, node2, data) = edge
        for id1 in id_map[node1]:
            for id2 in id_map[node2]:
                if id1 == id2:
                    continue
                if H.has_edge(id1, id2) and H[id1][id2] != data:
                    sys.stderr.write("Warning: Interaction %s <-> %s has a duplicate with different edge data\n" % (id1, id2))
                H.add_edge(id1, id2, data)
    return H

def readDB(conn, min_pubs=0):
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
    if min_pubs == 0:
        cursor.execute("SELECT entrez_id1, entrez_id2, count(pubmed_id) FROM interactions LEFT JOIN publications ON interactions.int_id = publications.int_id GROUP BY interactions.int_id")
    else:
        cursor.execute("SELECT MAX(counted) FROM (SELECT COUNT(pubmed_id) AS counted FROM publications GROUP BY int_id) AS x")
        max_pubs = cursor.fetchone()[0]
        cursor.execute("SELECT entrez_id1, entrez_id2, count(pubmed_id) FROM interactions LEFT JOIN publications ON interactions.int_id = publications.int_id GROUP BY interactions.int_id HAVING count(pubmed_id) >= %s", [min_pubs])
    for row in cursor.fetchall():
        G.add_edge(int(row[0]),int(row[1]),publications=int(row[2]))


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

    # Remove self-loops
    for node in G.nodes():
        try:
            G.remove_edge(node,node)
        except:
            continue

    return G
    
def readTSV(filename, delimiter="\t"):
    """Reads a Graph from a TSV formatted file.
    The TSV format should be as follows:
        Node1ID    Node2ID

    Parameters
    ----------
    filename : A file location

    delimiter : The delimiter for the file (defulat="\\t")

    Returns
    -------
    G : A NetworkX Graph
    """
    
    G = nx.Graph()
    with open(filename,'r') as file:
        next(file) # skip header
        for line in file:
            (idA, idB) = line.strip().split(delimiter)
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
            
            idA = fields.pop(0)
            type = fields.pop(0)
            
            for idB in fields:
                if idA != idB:
                    G.add_edge(idA, idB, type=type)

    return G

def writeGML(G, filename):
    """Writes a given network to file in GML format.

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
