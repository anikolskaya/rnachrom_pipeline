import sys
import pathlib
import pandas as pd
import argparse
import pyranges as pr
import os

from io_.fileops import load_rdc, load_blacklist
from pathlib import Path

parser = argparse.ArgumentParser(description='Removes contacts in which DNA parts overlap with ENCODE blacklisted regions')
parser.add_argument('path_rdc', type=str, help='Path to the rna-dna contacts file')
parser.add_argument('genome', help='Genome abbreviation (ENCODE blacklisted regions are available for hg19, hg38, mm9, mm10, dm3, ce10)')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')
args = parser.parse_args()

def remove_blacklisted_regions(path_rdc, genome):

    """
   Removes contacts in which DNA parts overlap with ENCODE blacklisted regions
    

    Parameters
    ----------
    rdc_path: str
        Path to an RDC-like file with RNA-DNA contacts
    genome: str 
        Genome abbreviation

    Returns
    -------
    pd.DataFrame 
    
    """
    cntcs = pr.PyRanges(load_rdc(path_rdc, sort = False))
    blist = pr.PyRanges(load_blacklist(genome))
    blck_removed = cntcs.overlap(blist, invert = True)

    calculate_stats(cntcs.as_df(), blck_removed.as_df())

    blck_removed.as_df().to_csv(os.path.join(args.outdir, Path(sys.argv[1]).stem + '_blacklst.tsv'), sep = '\t',
                                                                                                  header=False, 
                                                                                                  index=False)
    
    return 

def calculate_stats(before_blckl, after_blckl):
    """
    Returns file with contact counts just before and after removing blacklisted regions

    Parameters
    ----------
    before_blckl: pd.DataFrame
   
    after_blckl: pd.DataFrame

    outdir: str
    
    """
    stats = {}
    stats['raw'] = before_blckl.shape[0]
    stats['after_blacklisting'] = after_blckl.shape[0]
    pd.DataFrame([stats]).to_csv(os.path.join(args.outdir, Path(sys.argv[1]).stem + '.blacklisting.stat.tsv'), sep = '\t', 
                                                                                                  index=False)
    return 
    

pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
remove_blacklisted_regions(args.path_rdc, args.genome)
