#!/bin/python

import genes

def readBIFile(fname, conn=None):
    """Reads a pathway file as retrieved from the Broad Institute
    (http://software.broadinstitute.org/gsea/msigdb/collections.jsp)

    Note: Though these files are labeled as .gmt, they do not follow 
    the format whatsoever: 
    http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29

    Parameters
    ----------
    fname : the path to the broad institute label file

    conn : if given, the Entrez IDs will be validated against the genes database (optional)

    Returns
    -------
    paths : a dictionary with pathways as keys and sets of gene Entrez IDs as values
    """

    paths = {}
    eid_cache = {}

    with open(fname,'r') as file:
        for line in file:
            fields = line.strip().split("\t")

            pathway = fields[0]
            url = fields[1]
            genes = set(fields[2:])

            if conn is None: 
                paths[pathway] = genes
            else:
                eids = set()
                for gene in genes:
                    if gene not in eid_cache:
                        eid_cache[gene] = genes.checkEID(gene, conn)
                    eids.update(eid_cache[gene])
                paths[pathway] = eids

    return eids

def mapPaths(eids, paths):
    """Maps a given set of genes onto a set of paths.

    Parameters
    ----------
    eids : an iterable of Entrez IDs

    paths : a dictionary of paths to sets of Entrez IDs

    Returns
    -------
    matchedPaths : a dictionary of paths to Entrez IDs
    """

    eids_set = set(eids)
    return {path : eids_set & paths[path] for path in paths}
