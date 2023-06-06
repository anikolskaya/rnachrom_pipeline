import pyranges as pr
import pandas as pd
import ssl

from .schemas import rdc_pyranges_BED, rdc_dtypes, voted_BED, annot_pyranges_BED, BED3

def load_gtf_restricted(path):
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
        print("A GTF-like DataFrame is required")
        
    ann.gene_length = ann.lengths()
    ann = ann[ann.Feature == 'gene'][['gene_name', 'gene_type', 'gene_length']]
    
    return ann

def load_BED_annot(path):
    """
    Read intervals and metadata from a bed file.

    Parameters
    ----------
    path : str
        Path to a bed file

    Returns
    -------
    PyRanges

    """
    try:
            ann = pd.read_csv(path, 
                                          header = 0, 
                                          index_col=None,
                                          sep='\t',
                                          names = annot_pyranges_BED)
    except TypeError:
            print("A BED-like DataFrame is required")

    ann = pr.PyRanges(ann)
    ann.gene_length = ann.lengths()

    return ann[['gene_name', 'gene_type', 'gene_length']]
    

def load_rdc(path, header = None, sort = True, ncpus = 1):
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
    DataFrame

    """
    try:
        rdc = pd.read_csv(path,
                       header=header,
                       index_col=None,
                       sep='\t',
                       names=rdc_pyranges_BED,
                       dtype=rdc_dtypes)
    except ValueError:
        print("An RDC-like DataFrame is required")

    return rdc


def load_blacklist(genome):
    genomes = {'hg19': 'human', 'hg38': 'human',
                        'mm9': 'mouse', 'mm10': 'mouse'}
    ssl._create_default_https_context = ssl._create_unverified_context
    if genome != 'hg19':
        return pd.read_csv("https://mitra.stanford.edu/kundaje/akundaje/release/blacklists/{0}-{1}/{0}.blacklist.bed.gz".format(genome, genomes[genome]), sep = '\t', header = None,
	usecols = [0,1,2], names = BED3)
    else:
        return pd.read_csv("http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz", sep = '\t', header = None,
	usecols = [0,1,2], names = BED3)
    
