import pyranges as pr
import pandas as pd

from .schemas import rdc_pyranges_BED, rdc_dtypes, voted_BED

def load_gtf_restricted(path: str) -> pr.PyRanges:
    """
    Read intervals and metadata from a gtf file.

    Parameters
    ----------
    path : str
        Path to a gtf file

    Returns
    -------
    PyRanges

    """
    try:
        ann = pr.read_gtf(path)
    except TypeError:
        raise ("A GTF-like DataFrame is required")
        
    ann.gene_length = ann.lengths()
    ann = ann[ann.Feature == 'gene'][['gene_name', 'gene_type', 'gene_length']]
    
    return ann
    
    

def load_rdc(path: str, sort = True, ncpus = 1) -> pr.PyRanges:
    """
    Read intervals and metadata from a RNA-DNA contacts file.

    Parameters
    ----------
    path : str
        Path to a rdc file
        
    ncpus: int
        How many cpus to use. 
        Can at most use 1 per chromosome or chromosome/strand tuple. 
        Will only lead to speedups on large datasets.

    Returns
    -------
    PyRanges

    """
    try:
        rdc = pd.read_csv(path,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=rdc_pyranges_BED,
                       dtype=rdc_dtypes)
    except ValueError:
        raise ("An RDC-like DataFrame is required")
        
    rdc.drop(['rna_cigar', 'dna_cigar'], axis=1, inplace=True)
    rdc = pr.PyRanges(rdc)
    
    if sort:
        rdc = rdc.sort(nb_cpu = ncpus)
    
    return rdc