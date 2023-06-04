import pyranges as pr
import pandas as pd
import sys
import argparse
import os
from pathlib import Path

from io_ import *
from stats_voting import *
from io_.schemas import voted_BED, rdc_pyranges_BED
from io_.fileops import load_BED_annot
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Annotate RNA-parts of contacts, perform a voting procedure and deplete contacts derived from ribosomal RNAs')
parser.add_argument('rdc', type=str, help='Path to rna-dna contacts file')
parser.add_argument('gene_annot', type=str, help='Path to gene annotation')
parser.add_argument('--annot_file_format', nargs='?', default='GTP', type=str, help='Gene annotation file format: GTF or BED (default: GTP)')
parser.add_argument('--outdir', nargs='?', type=str, default = './results', help='A folder in which to store results')
#parser.add_argument('--ncpus', nargs='?', type=int, default=1, help='How many cpus to use')
args = parser.parse_args()


def run_annotation_and_voting(
                              gene_annot_path, 
                              contacts_path,
                              annot_format = 'GTF',
                              outdir = './results',
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
    
    ncpus: int
        How many cpus to use. 
        Can at most use 1 per chromosome or chromosome/strand tuple. 
        Will only lead to speedups on large datasets.


    Returns
    -------
    pd.DataFrame 
    
    """

    d = {}
    empty = []
    
    print(annot_format)
    if annot_format == 'GTF':
        gene_ann = load_gtf_restricted(gene_annot_path)
    elif annot_format == 'BED':
        gene_ann = load_BED_annot(gene_annot_path)
    else:
        raise Exception('Gene annotation format is not supported')
        
    cnts = load_rdc(contacts_path, header = None,
                    ncpus=ncpus)
    cnts.drop(['rna_cigar', 'dna_cigar'], axis=1, inplace=True)
    cnts = pr.PyRanges(cnts)
    
    d['selected_annot'], d['complement_annot'], d['no_annot'] = annotate_rdc(cnts, gene_ann, cpus = ncpus)
    
    d = dict(map(lambda item: (item[0], item[1].drop(like="annot$").as_df()), d.items()))
    
    d['voted_annot'] = vote(d['selected_annot'])
    print('Voting is completed..')
    
    d['voted_noribo'] = remove_ribo(d['voted_annot'])
    print('Ribosomal contacts are removed..')

    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    for name, dataframe in d.items():
        if not dataframe.empty:
            dataframe.iloc[:, :11].to_csv(os.path.join(outdir, Path(args.rdc).stem + '_' + name + '.tab'), sep = '\t',
                                                                                                                  header=voted_BED, 
                                                                                                                  index=False)
        else:
            dataframe.to_csv(os.path.join(outdir, Path(args.rdc).stem + '_' + name + '.tab'))
            empty.append(name)
    
    for empty_df in empty:
        d.pop(empty_df, None)

    calculate_stats(d)
    print('Stats are calculated..')

    return 

    

    
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

run_annotation_and_voting(args.gene_annot, args.rdc, args.annot_file_format, args.outdir)
