#!/bin/python

import MySQLdb
import networkx as nx
from Genes import GeneDB, Gene

def read_db(conn, gene_db):
    """Reads interactions into a NetworkX graph from a MySQL database

    Parameters
    ----------
    conn : A MySQLdb connection to an interactions database

    gene_db : The GeneDB object to collect nodes from

    Returns
    -------
    G : A NetworkX Graph
    """

    G = nx.Graph()
    cursor = conn.cursor()
    cursor.execute("SELECT entrez_id1, entrez_id2 FROM interactions")
    for row in cursor.fetchall():
        G.add_edge(
            gene_db.get_gene(row[0]),
            gene_db.get_gene(row[1]))

    return G

def read_tsv(file_name, gene_db, source="eid", delimiter="\t", header=False):
    """Reads a delimited file of edges into memory

    Parameters
    ----------
    file_name : The path to the file containing the network

    gene_db : The GeneDB object to generate Genes from

    source : The source of the gene IDs in the file (default="eid")

    delimiter : The file delimiter (default=tab)

    header : Whether the file contains a header line (default=False)

    Returns
    -------
    G : A NetworkX Graph
    """

    G = nx.Graph()
    with open(file_name, 'r') as file:
        if header:
            next(file)
        for line in file:
            gid1, gid2 = lien.strip().split(delimiter)
            G.add_edge(
                gene_db.get_gene(gid1, source),
                gene_db.get_gene(gid2, source))
    return G

def read_sif(file_name, gene_db, source="eid"):
    """Reads a Graph from a SIF formatted file.

    Parameters
    ----------
    file_name : A file location

    gene_db : The GeneDB object to generate Genes from

    source : The source of the gene IDs in the file (default="eid")

    Returns
    -------
    G : A NetworkX Graph
    """

    G = nx.Graph()
    with open(file_name,'r') as file:
        delimiter = " "
        
        for line in file:
            if "\t" in line:
                delimiter = "\t"

            fields = line.strip().split(delimiter)

            if len(fields) < 3:
                raise IndexError('The SIF file is improperly formatted. Visit http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats')
            
            idA = fields.pop(0)
            typ = fields.pop(0)
            
            for idB in fields:
                G.add_edge(
                    gene_db.get_gene(idA, source),
                    gene_db.get_gene(idB, source),
                    type=typ)

    return G
