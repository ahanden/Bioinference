#!/bin/python

from fisher import pvalue
import genes

def readDB(conn):
    """Reads annotations in from the genes database.

    Parameters
    ----------
    conn : A MySQLdb connection to the genes database.

    Returns
    -------
    paths : A dictionary of annotations/pathways to sets of genes

    Notes
    -----
    This method puts all annotations into the same key list - not separating
    by database/source. I may want to change this in the future.
    """

    paths = {}

    cursor = conn.cursor()
    cursor.execute("SELECT db, annotation, entrez_id FROM annotations")
    for row in cursor.fetchall():
        db, ann, eid = row
        label = "%s: %s" % (db, ann)
        if label not in paths:
            paths[label]  = set()
        paths[label].add(eid)

    return paths


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

def monteCarloEnrichment(eids, paths, iters, func, *args):
    """Computes pathway enrichment using Monte Carlo estimation.

    Parameters
    ----------
    eids : the genes to enrich

    paths : a dictionary of paths to sets of Entrez IDs

    iters : the number of iterations to perform for selection

    func : a function to select genes - it must return an
           iterable of genes.

    *args : arguments to pass to the function
    Returns
    -------
    results : a dictionary with the following structure:
              { pathway : {
                  pval,
                  corr_pval,
                  genes
              }
    """

    eids = set(eids)

    resuts = {}
    for path, pop in paths.iteritems():
        results[path] = {'genes' : eids * pop, 'pval' : 0}

    for i in xrange(iters):
        sample = set(func(*args))
        for path, pop in paths.iteritems():
            if len(sample & pop) >= len(results[path]['genes']):
                results[path]['pval'] += 1

    for path, vals in results.iteritems():
        vals['pval'] /= float(iters)

    fdrcorrection(results)

    return results

def fisherEnrichment(eids, paths, N=20000):
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
       
    fdrcorrection(results)

    return results

def fdrcorrection(results):
    """Applies Benjamini-Hochberg correction to a results set from
    one of this script's enrichment functions.

    This method alters the input directly and does not return anything.

    Parameters
    ----------
    results : a dictionary with pathways as keys and a dictionary containing
              a p value as values.
    """
    sorted_paths = sorted(results.keys(), key=lambda x: results[x]['pval'])
    sorted_pvals = sorted([results[path]['pval'] for path in results])

    nobs = float(len(sorted_pvals))
    ecdffactor = [x/nobs for x in xrange(1, int(nobs+1))]

    pvals_corrected_raw = [sorted_pvals[i] / ecdffactor[i] for i in xrange(len(sorted_pvals))]

    pvals_corrected =  [pvals_corrected_raw[-1]]
    pvals_corrected += [min(pvals_corrected_raw[i], pvals_corrected_raw[i+1]) for i in xrange(len(pvals_corrected_raw)-2,-1,-1)] 
    pvals_corrected =  pvals_corrected[::-1]
    pvals_corrected =  [p if p < 1 else 1 for p in pvals_corrected]

    for i in xrange(len(sorted_paths)):
        results[sorted_paths[i]]['corr_pval'] = pvals_corrected[i]

