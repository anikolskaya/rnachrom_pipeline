import re
import sys
import pathlib
import pandas as pd
import numpy as np
import argparse
import os

from io_.fileops import load_rdc
from pathlib import Path

parser = argparse.ArgumentParser(description='Parse CIGAR field and process spliced contacts')
parser.add_argument('path_rdc_primary', type=str, help='Path to the primary rna-dna contacts file')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')
args = parser.parse_args()

def run_cigar_processing(rdc_path, outdir = './results'):
    """
    Parsing CIGAR field and spliced contacts processing: remove complex splicing cases ('I', 'D', > 1 'N' in CIGAR field), retain longer part of spliced contact ('N' = 1)
    

    Parameters
    ----------
    rdc_path: str
        Path to an RDC-like file with RNA-DNA contacts

    Returns
    -------
    pd.DataFrame 
    
    """
    rdc = load_rdc(rdc_path, header = 0)
    rdc = splicing_cases(rdc)
    
    if not (rdc['N_cnt'] == 0).all() & (rdc['ID_cnt'] == 0).all():
        rdc = process_simple_splicing(rdc)
    
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    rdc.drop(['N_cnt', 'ID_cnt'], axis=1).to_csv(os.path.join(outdir, Path(sys.argv[1]).stem + '_CIGAR_processed.tsv'), sep = '\t',
                                                                                                  header=False, 
                                                                                                  index=False)
    return


def splicing_cases(df):
    df['N_cnt'] = df.rna_cigar.str.count('N')
    df['ID_cnt'] = df.rna_cigar.str.count(r'[ID]')
    return df


def process_simple_splicing(df):
    df_ss = df[(df.N_cnt == 1) & (df.ID_cnt == 0)]
    df_ss[['N1','N2']] = df_ss['rna_cigar'].str.split(r'[NM]', expand=True).iloc[:,[0,2]]
    df_ss = df_ss.astype({'N1':'int','N2':'int'})
    
    df_ss['rna_end'] = np.select([df_ss.N1 >= df_ss.N2], [df_ss.rna_bgn + df_ss.N1], default = df_ss.rna_end)
    df_ss['rna_bgn'] = np.select([df_ss.N1 < df_ss.N2], [df_ss.rna_end - df_ss.N2], default = df_ss.rna_bgn)
    df_ss.drop(['N1', 'N2'], axis = 1, inplace = True)
    
    return pd.concat([df[(df.N_cnt == 0) & (df.ID_cnt == 0)], df_ss])

run_cigar_processing(args.path_rdc_primary, args.outdir)
