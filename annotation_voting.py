import pyranges as pr
import pandas as pd
import numpy as np
import argparse
import warnings

from stats_voting import *
from io_.schemas import *
from io_.fileops import load_rdc, load_BED, make_pyranges_format, save_file_custom

warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Annotates RNA-parts of contacts, applies a voting procedure and removes contacts derived from ribosomal DNA')
parser.add_argument('rdc', type=str, help='Path to RNA-DNA contacts file')
parser.add_argument('annot', type=str, help='Path to the gene annotation')
parser.add_argument('--annot_format', nargs='?', default='GTF', type=str, choices = ['GTF', 'BED'], help='Gene annotation file format: GTF or BED (default: GTF)')
parser.add_argument('--no_ribo', action="store_true", help="Where to exclude contacts of rRNA")
parser.add_argument('--cpus', nargs='?', type=int, default = 1)
parser.add_argument('--outdir', nargs='?', type=str, default = './results', help='A folder to save results')

args = parser.parse_args()


def run_annotation_and_voting(gene_annot_path, contacts_path, annot_format, outdir, ncpus = 1):
    """
    Annotates RNA-parts of contacts, applies a voting procedure and removes
    contacts derived from ribosomal DNA.


    Parameters
    ----------
    gene_annot_path: str
        Path to the file with RNA annotation

    contacts_path: str
        Path to the file with RNA-DNA contacts

    annot_format: str
        GTF or BED

    outdir: str


    Returns
    -------
    pd.DataFrame

    """

    d = {}
    save_names = {'singletons': singleton_suffix,
                  'selected_annot': selected_ann_suffix,
                  'voted': voted_suffix,
                  'complement_annot': complement_ann_suffix,
                  'voted_noribo': voted_nr_suffix
                  }
    empty = []

    if annot_format == 'GTF':
        gene_ann = pr.read_gtf(gene_annot_path)

    elif annot_format == 'BED':
        gene_ann = load_BED(gene_annot_path)
        gene_ann.rename(columns={"name": "gene_name"}, inplace = True)
        gene_ann = make_pyranges_format(gene_ann, strand = True, save_old_names = False)
    else:
	raise Exception('Gene annotation format is not supported')

    gene_ann.gene_length = gene_ann.lengths()
    gene_ann = gene_ann[['gene_name', 'gene_type', 'gene_length', 'source']]

    cnts = load_rdc(contacts_path, header = 0, names = rdc_BED)
    cnts.drop(['rna_cigar', 'dna_cigar'], axis=1, inplace=True)
    cnts, old_cols = make_pyranges_format(cnts, strand = True)

    d['voted'], d['complement_annot'] = annotate_rdc(cnts, gene_ann, cpus = ncpus)
    del gene_ann

    d = dict(map(lambda item: (item[0], item[1].drop(like="annot$").as_df()), d.items()))

    d['voted'] = vote(d['voted'])

    cnts = cnts.as_df()
    cnts.colnames = old_cols

    d['singletons'] = cnts[~cnts['id'].isin(d['voted']['id'])]

    if args.no_ribo:
        d['voted_noribo'] = remove_ribo(d['voted'])

    for name, dfr in d.items():
        if not dfr.empty:
            if name == 'singletons':
                save_file_custom(dfr, outdir, args.rdc, save_names[name], hdr=voted_BED[:-2])

            else:
                save_file_custom(dfr.iloc[:,:12], outdir, args.rdc, save_names[name], hdr=voted_BED)

        else:
            save_file_custom(dfr, outdir, args.rdc, save_names[name], hdr=voted_BED)
            empty.append(name)

    del d['singletons']
    for empty_df in empty:
        d.pop(empty_df, None)

    calculate_stats_reduced(d, args.rdc, args.outdir)

    return

def annotate_rdc(contacts, annot, cpus):
    """
    Creates annotated genomic intervals from a RNA-DNA contacts file.

    Parameters
    ----------
    contacts: PyRanges
        RDC-like file converted to PyRanges

    annot: PyRanges
        A GTF- or BED-like file converted to PyRanges

    cpus: int


    Returns
    -------
    Three PyRanges:
         1. Intervals annotated by genes on selected strand
         2. Intervals annotated by genes on the complementing strand
         3. Intervals that did not overlap with the annotation

    """
    all_annot = contacts.join(annot, how = 'left', strandedness = False, suffix = '_annot', nb_cpu = cpus)
    strand_match = all_annot.Strand == all_annot.Strand_annot
    gene_found = all_annot.gene_name != '-1'

    return all_annot[(strand_match & gene_found)], all_annot[(~strand_match & gene_found)]

def vote(an_contacts):
    """
    Applies voting procedure: in the case of unambiguous annotation
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

    if ('source' in an_contacts.columns) and any('gencode' in s for s in an_contacts['source'].unique()):
        print('genecode has a priority')
        an_contacts['gencode'] = np.where(an_contacts['source'].str.contains('gencode'), 'gencode', 'other')
        return an_contacts.sort_values(['density', 'gencode'], ascending=(False, True)).drop_duplicates('id', keep='first').drop('gencode', axis=1)

    else:
	return an_contacts.sort_values('density', ascending=False).drop_duplicates('id', keep='first')

def remove_ribo(contacts, ribo_list = ['rRNA', 'rRNA_pseudogene', 'rRNA_RepM']):
    """
    Removes contacts mapping to ribosomal DNA


    Parameters
    ----------
    contacts: pd.DataFrame
        RDC-like file


    Returns
    -------
    pd.DataFrame

    """
    return contacts[~contacts.gene_type.isin(ribo_list)]

run_annotation_and_voting(gene_annot_path = args.annot,
                          contacts_path = args.rdc,
                          annot_format = args.annot_format,
                          outdir = args.outdir,
                          ncpus = args.cpus)
