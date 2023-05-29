import pyranges as pr
import pandas as pd
import sys
import argparse

from io_ import *
from stats_voting import *
from io_.schemas import voted_BED

parser = argparse.ArgumentParser(description='Annotate RNA-parts of contacts, perform a voting procedure and deplete contacts derived from ribosomal RNAs')
parser.add_argument('rdc', type=str, help='Path to rna-dna contacts file')
parser.add_argument('gene_annot', type=str, help='Path to gene annotation')
parser.add_argument('keep_strand', type=bool, help='Whether to overlap intervals on the same strand (True) or the opposite one (False)')
#parser.add_argument('ncpus',type=int, default=1, help='How many cpus to use')
args = parser.parse_args()


def run_annotation_and_voting(
                              gene_annot_path, 
                              contacts_path,
                              keep_strand = True,
                              ncpus = 1):
    """
    Annotate RNA-parts of contacts, perform a voting procedure and deplete
    contacts derived from ribosomal RNAs.
    

    Parameters
    ----------
    gene_annot_path: str
        Path to a GTF-like file with genomic annotation
    
    contacts_path: str
        Path to an RDC-like file with RNA-DNA contacts
        
    keep_strand: bool
        Whether to overlap intervals on the same strand (True) or the opposite one (False).
    
    ncpus: int
        How many cpus to use. 
        Can at most use 1 per chromosome or chromosome/strand tuple. 
        Will only lead to speedups on large datasets.


    Returns
    -------
    pd.DataFrame 
    
    """

    d = {}
    gene_ann = load_gtf_restricted(gene_annot_path)
    cnts = load_rdc(contacts_path, 
                    ncpus=ncpus)
    
    if keep_strand:
        d['selected_annot'], d['complement_annot'], d['no_annot'] = annotate_rdc(cnts, gene_ann, cpus = ncpus)
    else:
        d['complement_annot'], d['selected_annot'], d['no_annot'] = annotate_rdc(cnts, gene_ann, cpus = ncpus)
    
    d = dict(map(lambda item: (item[0], item[1].drop(like="annot$").as_df()), d.items()))
    
    d['voted_annot'] = vote(d['selected_annot'])
    
    d['voted_noribo'] = remove_ribo(d['voted_annot'])
    d['voted_noribo'].drop(['gene_length', 'count', 'density'], axis = 1, inplace = True)

    calculate_stats(d)
    
    d['voted_noribo'].columns = voted_BED
    
    return d

    

    
def annotate_rdc(contacts, annot, cpus = 1):
    """
    Create annotated genomic intervals from a RNA-DNA contacts file.

    Parameters
    ----------
    contacts: PyRanges
        RDC-like file converted to PyRanges
        
    annot: PyRanges
        A GTF-like file converted to PyRanges
        
    ncpus: int
        How many cpus to use. 
        Can at most use 1 per chromosome or chromosome/strand tuple. 
        Will only lead to speedups on large datasets.

    Returns
    -------
    Three PyRanges: 
         1. Intervals annotated by genes on selected strand
         2. Intervals annotated by genes on the complementing strand
         3. Intervals that did not overlap with the annotation

    """
    all_annot = contacts.join(annot, how = 'left', strandedness = False, 
                                    suffix = '_annot', nb_cpu = cpus)
    strand_match = all_annot.Strand == all_annot.Strand_annot
    gene_found = all_annot.gene_name != '-1'
    
    return all_annot[(strand_match & gene_found)], all_annot[(~strand_match & gene_found)], all_annot[~gene_found]



def vote(an_contacts):
    """
    Perform so-called voting procedure: in the case of unambiguous annotation 
    the preference is given to the gene with more dense coverage of RNA-parts
    

    Parameters
    ----------
    contacts: pd.DataFrame
        RDC-like file


    Returns
    -------
    pd.DataFrame 
    
    """
    an_contacts['count'] = an_contacts['gene_name'].map(an_contacts['gene_name'].value_counts())
    an_contacts['density'] = (an_contacts['count'] / 
                                 an_contacts['gene_length']) * 100000
    
    return an_contacts.sort_values('density', ascending=False).drop_duplicates('id', keep='first')



def remove_ribo(contacts,
                ribo_list = ['rRNA', 'rRNA_pseudogene', 'rRNA_RepM']):
    """
    Remove contacts belonging to ribosomal gene types
    

    Parameters
    ----------
    contacts: pd.DataFrame
        RDC-like file


    Returns
    -------
    pd.DataFrame 
    
    """
    return contacts[~contacts.gene_type.isin(ribo_list)]

run_annotation_and_voting(args.gene_annot, args.rdc, args.keep_strand)
