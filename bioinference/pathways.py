#!/bin/python

from fisher import pvalue
import genes

def readBIFile(fname, conn=None):
    """Reads a pathway file as retrieved from the Broad Institute
    (http://software.broadinstitute.org/gsea/msigdb/collections.jsp)

    Parameters
    ----------
    fname : the path to the broad institute label file

    conn : if given, the Entrez IDs will be validated against the genes database (optional)

    Returns
    -------
    paths : a dictionary with pathways as keys and sets of gene Entrez IDs as values
    
    Notes
    -----
    Though these files are labeled as .gmt, they do not follow the format whatsoever: 
    http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
    """

    paths = {}
    eid_cache = {}

    with open(fname,'r') as file:
        for line in file:
            fields = line.strip().split("\t")

            pathway = fields[0]
            url = fields[1]
            genes = set([int(x) for x in fields[2:]])

            if conn is None: 
                paths[pathway] = genes
            else:
                eids = set()
                for gene in genes:
                    if gene not in eid_cache:
                        eid_cache[gene] = genes.checkEID(gene, conn)
                    eids.update(eid_cache[gene])
                paths[pathway] = eids

    return paths

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

    eids = set(eids)
    return {path : eids & paths[path] for path in paths}

def enrichPathways(eids, paths, N=20000):
    """Computes pathway enrichment using Fisher's Exact Test
    with FDR correction using the Benjamini-Hochberg procedure.

    Parameters
    ----------
    eids : the genes to enrich

    paths : a dictionary of paths to sets of Entrez IDs

    N : the population size (default=20000)

    Returns
    -------
    results : a dictionary with the following structure:
              { pathway : {
                  pval,
                  corr_pval,
                  genes
              }
    """
    sample = set(eids)

    results = {}

    for path, pop in paths.iteritems():
        tp = len(sample & pop)
        fp = len(sample) - tp
        fn = len(pop) - tp
        tn = N - tp - fp - fn
        
        p = pvalue(tp, fp, fn, tn).right_tail

        results[path] = {'pval': p, 'genes': sample & pop}
        
    sorted_paths = sorted(results.keys(), key=lambda x: results[x]['pval'])
    sorted_pvals = sorted([results[path]['pval'] for path in results])

    nobs = float(len(sorted_pvals))
    ecdffactor = [x/nobs for x in xrange(1, int(nobs+1))]

    pvals_corrected_raw = [sorted_pvals[i] / ecdffactor[i] for i in xrange(len(sorted_pvals))]

    pvals_corrected = [pvals_corrected_raw[-1]]
    pvals_corrected += [min(pvals_corrected_raw[i], pvals_corrected_raw[i+1]) for i in xrange(len(pvals_corrected_raw)-2,-1,-1)] 
    pvals_corrected = pvals_corrected[::-1]
    pvals_corrected = [p if p < 1 else 1 for p in pvals_corrected]

    for i in xrange(len(sorted_paths)):
        results[sorted_paths[i]]['corr_pval'] = pvals_corrected[i]

    return results
