#!/usr/bin/python

from fisher import pvalue
from math import log
from Genes import GeneDB

def read_ts(file_name, gene_db, species_names=None):
    """Reads the Target Scan database into memory as a dictionary.

    Parameters
    ----------
    fle_name : The path to the TargetScan database dump

    gene_db : A GeneDB object to generate Genes from

    species_names: A list of the proper names of species to extract miRNAs
                  for. If given, this method will only return miRs for the given 
                  species. (Optional)

                  The following species names are supported:
                    Homo sapiens
                    Mus musculus
                    Rattus norvegicus
                    Monodelphis domestica
                    Xenopus tropicalis
                    Gallus gallus
                    Macaca mulatta
                    Pan troglodytes
                    Canis familiaris
                    Bos tarus
    Returns
    -------
    targets: A dictionary with miRs as keys and sets of genes as values.


    Note
    ----
    This method is currently only set up for human and mouse
    """

    targets = {}

    prefix_map = {
        9606:  "hsa",
        10090: "mmu",
        10116: "rno",
        13616: "mdo",
        8364:  "xtr",
        9031:  "gga",
        9544:  "mml",
        9598:  "ptr",
        9615:  "cfa",
        9913:  "bta"
    }

    species_map = {
        9606: 'Homo sapiens',
        10090: 'Mus musculus',
        10116: 'Rattus norvegicus',
        13616: 'Monodelphis domestica',
        8364: 'Xenopus tropicalis',
        9031: 'Gallus gallus',
        9544: 'Macaca mulatta',
        9598: 'Pan troglodytes',
        9615: 'Canis familiaris',
        9913: 'Bos tarus'
    }

    if species_names:
        species_names = set(species_names)

    with open(file_name, 'r') as file:
        # Skip the header
        next(file)

        for line in file:
            fields = line.strip().split("\t")

            mirs          = fields[0][4:].split("/")
            target_symbol = fields[2]
            species_code  = int(fields[4])

            if not species_names or species_map[species_code] in species_names:

                for mir in mirs:

                    mir = prefix_map[species_code]+"-miR-"+mir
                    if '.' in mir:
                        mir = mir[:mir.index('.')]

                    if mir not in targets:
                        targets[mir] = set()
                        
                    targets[mir].add(gene_db.get_gene(fields[2], "symbol"))

    return targets


def read_mtb(file_name, gene_db, species_names=None):
    """Reads the miRTarBase database into memory as a dictionary.

    Parameters
    ----------
    file_name : The path to the miRTarBase database dump

    gene_db : A GeneDB object to generate Genes from

    species_names: A list of the proper names of species to extract miRNAs
                  for. If given, this method will only return miRs for the given 
                  species. (Optional)

    Returns
    -------
    A dictionary with miRs as keys and sets of genes as values.
    """

    targets = {}

    if species_names:
        species_names = set(species_names)

    with open(file_name, 'r') as file:
        for line in file:
            fields = line.strip().split(",")

            mir     = fields[1]
            species = fields[2]
            target  = fields[4]

            if not species_names or species in species_names:

                if mir not in targets:
                    targets[mir] = set()

                targets[mir].update(gene_db.get_gene(target))

    return targets

def spanning_scores(mirs, graph, interactome, clusters):
    """Computes spanning scores for all miRNA.

    Parameters
    ----------
    mirs : A dictionary of miRs to their Gene targets

    graph : The graph over which to compute the spanning scores

    interactome; The consolidate interactome

    clusters : A dictionary of clusters for the graph

    Returns
    -------
    A dictionary of miRs to their corresponding spanning scores
    """

    stats = {}
    N = len(interactome)  # population size
    M = len(graph) # sample size

    clusters = [set(genes) for genes in clusters.values()]

    data = {}
    for mir, targets in mirs.iteritems():
        targets = set(targets)

        pool = targets & set(interactome) # get target pool
        n    = len(pool)                  # successes in population
        x    = len(pool & set(graph))     # successes in sample

        p_val = fishers(N,n,M,x)
        neg_p = -log(p_val)
        
        clust_count = sum([1 for clu in clusters if len(clu) > 3 and set(targets) & clu])

        data[miR] = {
            'hits'        : x,
            'neg_p'       : neg_p,
            'clust_count' : clust_count,
            'p_val'       : pval
        }

    max_o = float(max([0] + [data[miR]['hits']        for miR in data])) # best achieved overlap
    max_c = float(max([0] + [data[miR]['clust_count'] for miR in data])) # best achieved cluster count
    max_p = float(max([0] + [data[miR]['neg_p']       for miR in data])) # best achieved -log(pval)

    if max_o == 0:
        raise Exception("There are no miRNA targets in the network")
    if max_c == 0:
        raise Exception("There are no miRNA targets in any of the clusters")

    scores = {}
    for miR, stats in data.iteritems():

        first  = stats['negp'] / max_p
        second = stats['clust_count'] / max_c
        third  = stats['hits'] / max_o

        scores[miR] = first + second + third # sum score components

    return scores

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

