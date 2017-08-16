import MySQLdb
from fisher import pvalue
from math import log
import re
import sys
import genes


def loadTargetScan(conn, fname, species_names=None):
    """Reads the Target Scan database into memory as a dictionary.

    Parameters
    ----------
    conn : a MySQLdb connection to the genes database

    fname : The path to the TargetScan database dump

    species_names: A list of the proper names of species to extract miRNAs
                  for. If given, this method will only return miRs for the given 
                  species. (Optional)

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

    with open(fname, 'r') as file:
        next(file)

        for line in file:
            fields = line.strip().split("\t")

            mirs          = fields[0][4:].split("/")
            target_symbol = fields[2]
            species_code  = int(fields[4])

            if species_names and species_map[species_code] in species_names:

                for mir in mirs:

                    mir = prefix_map[species_code]+"-miR-"+mir
                    if '.' in mir:
                        mir = mir[:mir.index('.')]

                    if mir not in targets:
                        targets[mir] = set()
                        
                    targets[mir].update(genes.getEID(fields[2], conn))

    return targets


def loadMirTarBase(conn, fname, species_names=None):
    """Reads the miRTarBase database into memory as a dictionary.

    Parameters
    ----------
    conn : a MySQLdb connection to the genes database

    fname : The path to the miRTarBase database dump

    species_names: A list of the proper names of species to extract miRNAs
                  for. If given, this method will only return miRs for the given 
                  species. (Optional)

    Returns
    -------
    targets: A dictionary with miRs as keys and sets of genes as values.
    """

    targets = {}

    if species_names:
        species_names = set(species_names)

    with open(fname, 'r') as file:
        for line in file:
            fields = line.strip().split(",")

            mir     = fields[1]
            species = fields[2]
            target  = fields[4]

            if species_names and species in species_names:

                if mir not in targets:
                    targets[mir] = set()

                targets[mir].update(genes.checkEID(target, conn))

    return targets

def spanningScores(TS, NET, CI, clusters):
    """Computes spanning scores for all miRNA.

    Parameters
    ----------
    TS : the TargetScan dictionary

    NET : the network to compute scores over

    CI : the consolidated interactome

    clusters : a list of clusters for NET

    Returns
    -------
    """

    stats = {}
    N = len(CI)  # population size
    M = len(NET) # sample size

    data = {}
    for miR, targets in TS.iteritems():
        targpool = set(targets) & set(CI) # get target pool

        n = len(targpool)                # successes in population
        x = len(targpool & set(NET))    # successes in sample

        pval = fishers(N,n,M,x)
        negp = -log(pval)
        
        clust_count = sum([1 for cluster in clusters if len(cluster) >3 and set(targets) & set(cluster)])

        data[miR] = {
            'hits'        : x,
            'negp'        : negp,
            'clust_count' : clust_count,
            'pval'        : pval
        }

    max_o = float(max([0] + [data[miR]['hits']        for miR in data])) # best achieved overlap
    max_c = float(max([0] + [data[miR]['clust_count'] for miR in data])) # best achieved cluster count
    max_p = float(max([0] + [data[miR]['negp']        for miR in data])) # best achieved -log(pval)

    if max_o == 0:
        raise Exception("There are no miRNA targets in the network")
    if max_c == 0:
        raise Exception("There are no miRNA targets in any of the clusters")

    scores = {}
    for miR, stats in data.iteritems():

        first = stats['negp'] / max_p    # -log(pval) / -log(best_pval) (max = 1)
        second = stats['clust_count'] / max_c # cluster_count / best_cluster_count (max = 1)
        third = stats['hits'] / max_o # hit_count / best_hit_count (max = 1)

        scores[miR] = first + second + third    # sum score components (max score = 3)

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

