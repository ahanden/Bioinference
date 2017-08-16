#!/bin/python

from fisher import pvalue
from Genes import GeneDB

def read_bi_file(file_name, gene_db conn=None):
    """Reads a pathway file as retrieved from the Broad Institute
    (http://software.broadinstitute.org/gsea/msigdb/collections.jsp)

    Parameters
    ----------
    fname : The path to the broad institute label file

    gene_db : A GeneDB object to generate Genes with

    Returns
    -------
    A dictionary with pathways as keys and sets of gene Genes as values
    
    Notes
    -----
    Though these files are labeled as .gmt, they do not follow the format whatsoever: 
    http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
    """

    paths = {}
    eid_cache = {}

    with open(file_name,'r') as file:
        for line in file:
            fields  = line.strip().split("\t")

            pathway = fields[0]
            url     = fields[1]
            eids    = set([int(x) for x in fields[2:]])

            paths[pathway] = set([ gene_db.get_gene(eid) for eid in eids ])

    return paths


def fisher_enrichment(genes, paths, N=20000):
    """Computes pathway enrichment using Fisher's Exact Test
    with FDR correction using the Benjamini-Hochberg procedure.

    Parameters
    ----------
    genes : An interable of Gene objects to enrich

    paths : a dictionary of paths to Gene objects

    N : the population size (default=20000)

    Returns
    -------
    results : a dictionary with the following structure:
              { pathway : {
                  pval,
                  corr_pval,
                  genes,
                  size
              }
    """

    sample = set(genes)

    results = {}

    for path, pop in paths.iteritems():
        tp = len(sample & pop)
        fp = len(sample) - tp
        fn = len(pop) - tp
        tn = N - tp - fp - fn
        
        p = pvalue(tp, fp, fn, tn).right_tail

        results[path] = {
            'pval'  : p, 
            'genes' : sample & pop,
            'size'  : len(pop)}
       
    fdr_correction(results)

    return results

def fdr_correction(enrichments):
    """Applies Benjamini-Hochberg correction to enrichment enrichments.

    Parameters
    ----------
    enrichments : a dictionary with pathways as keys and a dictionary
                  containing a p-value as values.

    Returns
    -------
    This method does not return anything, but rather alters the enirchments
    dictionary in place. by adding a 'corr_pval' key to each value.
    """
    sorted_paths = sorted(enrichments.keys(), key=lambda x: enrichments[x]['pval'])
    sorted_pvals = sorted([enrichments[path]['pval'] for path in enrichments])

    nobs = float(len(sorted_pvals))
    ecdffactor = [x/nobs for x in xrange(1, int(nobs+1))]

    pvals_corrected_raw = [sorted_pvals[i] / ecdffactor[i] for i in xrange(len(sorted_pvals))]

    pvals_corrected =  [pvals_corrected_raw[-1]]
    pvals_corrected += [min(pvals_corrected_raw[i], pvals_corrected_raw[i+1]) for i in xrange(len(pvals_corrected_raw)-2,-1,-1)] 
    pvals_corrected =  pvals_corrected[::-1]
    pvals_corrected =  [p if p < 1 else 1 for p in pvals_corrected]

    for i in xrange(len(sorted_paths)):
        enrichments[sorted_paths[i]]['corr_pval'] = pvals_corrected[i]

