import sys
import argparse
import pandas as pd
import numpy as np
import os 
import pathlib
import concurrent.futures

from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

pd.options.mode.chained_assignment = None

from io_.schemas import chrsizes_format, normalized_BED
from io_.fileops import load_rdc


parser = argparse.ArgumentParser(description='Calculates inner- and inter-chromosomal scaling based on mRNA contacts density decay')
parser.add_argument('path_rdc_norm', type=str, help='Path to the rna-dna contacts file, background_normalized')
parser.add_argument('path_chrsizes', type=str, help='Path to the file with chromosome sizes')
parser.add_argument('--factor', type=float, nargs='?', default = 1.2, help='Scale factor')
parser.add_argument('--threads', type=int, nargs='?', default = 1)
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to save results')
args = parser.parse_args()

def run_all_scaling(cnts_all_path, chrsizes_path, b0 = 0, b_n = 100, factor = 1.2,
                threads = 1):

    cnts_all = load_rdc(cnts_all_path, sort = False, header = None, names = normalized_BED)

    chrsizes = pd.read_csv(chrsizes_path, sep = '\t', 
                           header = None, names = chrsizes_format)
    
    bins_ = make_geom_bins(length = max(chrsizes.len))
    
    #inner scaling, runs in parallel for each chromosome
    with concurrent.futures.ThreadPoolExecutor(threads) as pool:
        
        scaled_cis = pd.concat(pool.map(partial(run_inner_scaling, bins_), 
                                   [group for name, group in cnts_all.groupby('rna_chr')]))
    
    scaled_cis.to_csv(os.path.join(args.outdir, Path(cnts_all_path).stem + '_SCinnerCHR.tab'), index = False, sep = '\t')

    #inter scaling, runs in parallel for each chromosome
    with concurrent.futures.ThreadPoolExecutor(threads) as pool:
        
        scaled_trans = pd.concat(pool.map(partial(run_inter_scaling, chrsizes), 
                                   [group for name, group in cnts_all.groupby('rna_chr')]))
        
    scaled_trans.to_csv(os.path.join(args.outdir, Path(cnts_all_path).stem + '_SCinterCHR.tab'), index = False, sep = '\t')
        
    columns = ['bin', 'bin_size', 'bin_mean', 'n_PC_bin', 'dens_PC_bin', 'log10Dist_PC',
               'log10Dens_PC', 'n_PC_trans', 'dens_PC']
    
    pd.concat([scaled_cis, scaled_trans]).drop(columns = columns).to_csv(os.path.join(args.outdir, Path(cnts_all_path).stem + '_SC.tab'), index = False, sep = '\t')
    
    return 



def make_geom_bins(length, b0 = 0, b_n = 100, factor = 1.2):
    """
    Makes geometrically increasing bins

    Parameters
    ----------
    length : int
        Chromosome length
        
    b0 : the 1st term of a geometric sequence
    
    b_n : The first step of a geometric sequence
    
    factor: the common ratio

    Returns
    -------
    np.array
    
    """
    
    bins = [0]
    
    while b_n < length:
        bins.append(b_n)
        b_n = b0 + factor * b_n
        
    bins.append(length)
    
        
    return bins


def rna_dna_distance_binning(cnts, bins):
    
    """
    Assigns bin according to distance between RNA and DNA part of a contact

    Parameters
    ----------
    cnts : pd.DataFrame
        
    bins : list

    Returns
    -------
    pd.DataFrame

    """
    cnts = cnts.query("dna_chr == rna_chr")
    
    cnts['dist'] = (cnts['rna_bgn'] - cnts['dna_bgn']).abs()
    
    cnts['bin'] = pd.cut(cnts['dist'], bins = bins, include_lowest = True, 
                   precision = 2, duplicates = 'drop')

    cnts['bin_size'] = pd.to_numeric(cnts['bin'].apply(lambda x: (x.right - x.left)))
    cnts['bin_mean'] = pd.to_numeric(cnts['bin'].apply(lambda x: x.mid))
    
    return cnts


def inner_PC_density(cnts_b):
    
    """
    Calculates protein-coding RNA contacts density in each bin, based on RNA-DNA distance

    Parameters
    ----------
    cnts : pd.DataFrame


    Returns
    -------
    pd.DataFrame
        Contains computed values based on mRNA cis-contacts, assigned to each bin

    """
    cnts_b = cnts_b.query('gene_type == "protein_coding"')
    libS_pc = cnts_b['N2'].sum()
        
    cnts_b['n_PC_bin'] = cnts_b.groupby(by = ['bin'])['N2'].transform('sum')
    cnts_b['dens_PC_bin'] = (cnts_b['n_PC_bin'] / cnts_b['bin_size']) / libS_pc
        
    cnts_b['log10Dist_PC'] = np.log10(cnts_b['bin_mean'])
    cnts_b['log10Dens_PC'] = np.log10(cnts_b['dens_PC_bin'])
    
    
    bins_pc_stat = cnts_b.sort_values('bin').drop_duplicates('bin', keep='first')
        
    return bins_pc_stat[['bin', 'n_PC_bin', 'dens_PC_bin','log10Dist_PC', 'log10Dens_PC']]


def inter_PC_density(cnts, chrsizes):
    """
    Calculates protein-coding RNA contacts density on non-parental chromosomes

    Parameters
    ----------
    cnts : pd.DataFrame
    
    chrsizes: pd.DataFrame


    Returns
    -------
    pd.DataFrame
        Contains computed values based on mRNA trans-contacts, in chromosome-chromosome fashion

    """
    
    #select only mRNA trans-contacts
    cnts = cnts.query('dna_chr != rna_chr and gene_type == "protein_coding"')
    
    cnts['n_PC_trans'] = cnts.groupby(by = ['rna_chr', 'dna_chr'])["N2"].transform('sum')
    cnts = cnts.merge(chrsizes, left_on = 'dna_chr', right_on = 'chr')
    cnts['dens_PC'] = (cnts['n_PC_trans'] / cnts['len']) * 1000
    
    cnts.drop(columns=['chr'], inplace = True)
    
    bins_pc_stat = cnts.sort_values('dna_chr').drop_duplicates('dna_chr', keep='first')
    
    return bins_pc_stat[['dna_chr', 'n_PC_trans', 'dens_PC']].reset_index(drop=True)


def inner_scaling(cnts_cis, dens_pc_cis):
    """
    Calculates scaling for all RNA-DNA cis-contacts

    Parameters
    ----------
    cnts_cis : pd.DataFrame
    
    dens_pc_cis: pd.DataFrame
         Pre-calculated table from the function inner_PC_density()


    Returns
    -------
    pd.DataFrame
        
    """
    
    cnts_cis = cnts_cis.query('rna_chr == dna_chr')
    libS_cis = cnts_cis['N2'].sum()
    
    cnts_cis = cnts_cis.merge(dens_pc_cis, on = 'bin')
    
    cnts_cis['normSC_raw'] = (cnts_cis['N2'] / cnts_cis['dens_PC_bin']) 
    cnts_cis['normSC'] = cnts_cis['normSC_raw'] * (libS_cis / sum(cnts_cis['normSC_raw']))
    
    cnts_cis.drop(columns=['normSC_raw'], inplace = True)
    
    return cnts_cis


def inter_scaling(cnts_trans, dens_pc_trans):
    """
    Calculates scaling for all RNA-DNA trans-contacts

    Parameters
    ----------
    cnts_trans : pd.DataFrame
    
    dens_pc_trans: pd.DataFrame
        Pre-calculated table from the function inter_PC_density()


    Returns
    -------
    pd.DataFrame
        
    """
    
    cnts_trans = cnts_trans.query('rna_chr != dna_chr')
    libS_trans = cnts_trans['N2'].sum()
    
    cnts_trans = cnts_trans.merge(dens_pc_trans, on = 'dna_chr')
    
    cnts_trans['normSC_raw'] = (cnts_trans['N2'] / cnts_trans['dens_PC']) 
    cnts_trans['normSC'] = cnts_trans['normSC_raw'] * (libS_trans / sum(cnts_trans['normSC_raw']))
    
    cnts_trans.drop(columns=['normSC_raw'], inplace = True)
    
    return cnts_trans


def run_inner_scaling(bins, cnts_all):
    """
    Runs all commands to calculate inner scaling

    Parameters
    ----------
    bins : list
    
    cnts_all: pd.DataFrame

    Returns
    -------
    pd.DataFrame
        
    """
    cnts_inner_all = rna_dna_distance_binning(cnts_all, bins)
    #dataset with cis-contacts only, RNA-DNA distances are binned
    
    pc_values = inner_PC_density(cnts_inner_all)
    #scaling parameters
    
    pc_values.to_csv(os.path.join(args.outdir, Path(args.path_rdc_norm).stem + '_PCinnerCHR.tab'), index = False, sep = '\t')
    
    cnts_inner_scaled = inner_scaling(cnts_inner_all, pc_values)
    #scaling for all RNAs
    
    
    return cnts_inner_scaled


def run_inter_scaling(chrsizes, cnts_all):
    """
    Runs all commands to calculate inner scaling

    Parameters
    ----------
    bins : list
    
    cnts_all: pd.DataFrame

    Returns
    -------
    pd.DataFrame
        
    """

    pc_values = inter_PC_density(cnts_all, chrsizes)
    #scaling parameters
    
    pc_values.to_csv(os.path.join(args.outdir, Path(args.path_rdc_norm).stem + '_PCinterCHR.tab'), index = False, sep = '\t')
    
    cnts_inter_scaled = inter_scaling(cnts_all, pc_values)
    #scaling for all RNAs
    
    
    return cnts_inter_scaled


run_all_scaling(args.path_rdc_norm, args.path_chrsizes, b0 = 0, b_n = 100, factor = args.factor,
                threads = args.threads)
