import MySQLdb
from fisher import pvalue
from math import log
import re
import sys

def readTargetScanFile(conn, filename, species_code=9606):
    """Reads data from a Target Scan output file into memory.

    Parameters
    ----------
    conn : a MySQL connection to the genes database

    filename : A file location
    
    species_code : Which species to filter for (default=9606).
        Will not filter for species if the value is None.

    Returns
    -------
    TS : A dictionary of lists
        The dictionary keys are miRNAs, and values are lists of targets as Entrez IDs
    """
    
    TS = {}
    cursor = conn.cursor()
    p = re.compile('(ENSG\d+)\.?')

    with open(filename) as file:
        next(file) # Skip the header
        for line in file:
            fields = line.strip().split("\t")

            mirna         = fields[0]
            target_id     = fields[1]
            target_symbol = fields[2]
            transcript    = fields[3]
            species       = int(fields[4])
           
            # Skip along if this doesn't meet our filters
            if species_code is not None and species != species_code:
                continue

            eids = getEID(cursor, target_symbol)
            if not eids:
                m = p.match(target_id)
                eids = crossQuery(cursor, 'Ensembl', m.group(1)) 

            if eids:
                if mirna not in TS:
                    TS[mirna] = {}
                TS[mirna].update({eid: True for eid in eids})
    
    # Clean up formatting
    for mirna in TS:
        TS[mirna] = TS[mirna].keys()
    
    return TS
    
    
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

    max_o = float(max([data[miR]['hits']        for miR in data])) # best achieved overlap
    max_c = float(max([data[miR]['clust_count'] for miR in data])) # best achieved cluster count
    max_p = float(max([data[miR]['negp']        for miR in data])) # best achieved -log(pval)

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

