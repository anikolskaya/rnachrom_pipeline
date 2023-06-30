import argparse
import sys
import pandas as pd
import subprocess
import os 
import pathlib
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

from io_.fileops import load_rdc, load_BED
from io_.schemas import voted_BED, BED_Stereogene

parser = argparse.ArgumentParser(description='Calculate N2 background normalization')
parser.add_argument('path_rdc', type=str, help='Path to the rna-dna contacts file')
parser.add_argument('path_Stereogene', type=str, help='Path to the Stereogene smoother script')
parser.add_argument('path_config', type=str, help='Path to the config template file')
parser.add_argument('path_chrsizes', type=str, help='Path to the file with chromosome sizes')
parser.add_argument('--filter_top_n', type=int, nargs='?', default = 50, help='Drop the most contacting RNAs (rank, default: 50)')
parser.add_argument('--filter_tail_n', type=int, nargs='?', default = 1000, help='Drop the least contacting RNAs (rank, default: 1000)')
parser.add_argument('--outdir', type=str, nargs='?', default = './results', help='A folder in which to store the results')
parser.add_argument('--stat_folder', type=str, nargs='?', default = './stats', help='A folder in which to store stats')
args = parser.parse_args()


def run_normalization(rdc_path, stereogene_path, 
                      config_path, chrsizes_path,
                      top_n_drop, tail_n_drop, outdir, stats_outdir):
    
    """
    Runs all commands

    Parameters
    ----------
    rdc_path: str
        Path to the file with RNA-DNA contacts
        
    stereogene_path: str
        Path to the Stereogene Smoother
        
    config_path: str
        Path to the config template file
        
    chrsizes_path: str
        Path to the file with chromosome sizes
        
    top_n_drop: int
        Drop the most contacting RNAs (a rank)
        
    tail_n_drop: int
        Drop the least contacting RNAs (a rank)
        
    outdir: str
    
    """
    
    rdc = load_rdc(rdc_path, header = 0, 
                            sort = False, names = voted_BED)

    print('rdc is loaded')
    
    rdc_back = select_bg_rnas(rdc, filter_top_n = top_n_drop, 
                              filter_tail_n = tail_n_drop)

    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    rdc_back[['dna_chr', 'dna_bgn', 'dna_end']].to_csv(os.path.join(outdir, Path(rdc_path).stem + '_bg.bed'), sep = '\t', 
                                                                                                            index=False, header=False)
    
    mod_cfg = os.path.dirname(config_path) + '/cfg_mod.cfg'
    
    prepeare_cfg(cfg_path = config_path, 
                 chrsizes_path = chrsizes_path,
                 mod_cfg_path = mod_cfg)
    
    
    cmd = prepare_cmd(path_smoother = stereogene_path,
                  path_cfg = mod_cfg,
                  path_back = rdc_path,
                  outdir = outdir)

    print(cmd)
    
    run_smoother(cmd)
    
    bg_sm_path = os.path.abspath(os.path.join(outdir, Path(rdc_path).stem + '_bg_sm.bgr'))
    
    bg = backannot(path_to_bg = bg_sm_path, 
                   rdc = rdc)

    print("last step")
    
    calculate_N2_normalization(bg)

    back_stats(rdc, rdc_back).to_csv(os.path.join(stats_outdir, Path(rdc_path).stem + '_bg.tab'),  sep = '\t', 
                                                                                                            index=False)
    
    return




def select_bg_rnas(all_cnts, filter_top_n, filter_tail_n):
    """
    Prepares contacts for background model

    Parameters
    ----------
    rdc_path : pd.DataFrame
        Path to a rdc file
        
    filter_top_n: int
        Drop contacts from n the most high-contacting RNAs
        
    filter_tail_n:
        Drop contacts from n the most low-contacting RNAs

         
    Returns
    -------
    pd.Dataframe with selected background contacts
         
    """
    
    
    back_pc_rnas = all_cnts[all_cnts['gene_type'] == 'protein_coding']['gene_name'].value_counts()[filter_top_n:-filter_tail_n].index
    rdc_back = all_cnts.query("gene_name in @back_pc_rnas & rna_chr != dna_chr")
    
    return rdc_back


def back_stats(all_cnts,
               back_cnts):
    """
    Calculates statistics on background contacts

    Parameters
    ----------
    all_cnts: pd.DataFrame
        
    back_cnts: pd.DataFrame
        Table with background contacts
        
    Returns
    -------
    pd.DataFrame with statistics on background contacts
    
    """
    
    stats = {}
    stats['Total_cnts'] = all_cnts.shape[0]
    stats['Background_cnts'] = back_cnts.shape[0]
    stats['Background_cnts, %'] = round(back_cnts.shape[0] / all_cnts.shape[0] * 100, 3)
    
    return pd.DataFrame([stats])



def prepeare_cfg(cfg_path, chrsizes_path, mod_cfg_path):
    """
    Prepares Stereogene config file

    Parameters
    ----------
    cfg_path : str
        Path to template config file
        
    chrsizes_path: str
        Path to the file with chromosome sizes
        
    mod_cfg_path: str

    outdir: str
         Path to the output directory
    
    """
    
    with open(cfg_path, 'r') as raw_cfg:
        content = raw_cfg.read()
        
        for name, value in [('CHRSIZE', chrsizes_path), ('REPLACE', args.outdir)]:
            content = content.replace(name, value)

    with open(mod_cfg_path, 'w+') as mod_cfg:
        mod_cfg.write(content)
    
    return 


def prepare_cmd(path_smoother,
                path_cfg,
                path_back, outdir):
    
    """
    Prepares bash command to run Smoother

    Parameters
    ----------
    path_smoother : str
        Path to the Stereogene Smoother
        
    path_cfg: str
        Path to the modified config
        
    path_back: str
        Path to the background contacts file
         
    Returns
    -------
    List with a command
    
    """
    
    args_for_sg = [path_smoother,
        os.path.abspath(path_cfg),
        os.path.abspath(os.path.join(outdir, Path(path_back).stem + '_bg.bed'))]
    
    cmd = "{0} cfg={1} {2}".format(args_for_sg[0], args_for_sg[1], args_for_sg[2])
    return [cmd]

def run_smoother(cmd):
    """
    Runs Smoother

    Parameters
    ----------
    cmd : list
        A command to run
        
    Returns
    -------
    Saves the result in folder specified in config file
    
    """
    
    sb = subprocess.Popen(cmd, shell=True, executable='/bin/bash', text=True)
    
    sb.wait()
    
    if sb.returncode != 0:
        print('Non-zero exit status')
        return sb.communicate()[0]
    
    return 

def backannot(path_to_bg,
              rdc):
    """
    Annotate contacts with background weights

    Parameters
    ----------
    path_to_bg: str
        Path to the file with smoothed background counts
        
    rdc: pd.DataFrame
        Contacts to annotate
        
    Returns
    -------
    pd.DataFrame with weights which are written in bg_sm column
    
    """
    sm_bg = load_BED(path_to_bg, names = BED_Stereogene)
    sm_bg = pr.PyRanges(sm_bg)
    
    #place dna parts before rna parts 
    dna_first = rdc.iloc[:, np.r_[5:11, 0:5]]
    dna_first = pr.PyRanges(dna_first.rename(columns={"dna_chr": "Chromosome",
                              "dna_bgn": "Start",
                              "dna_end": "End"})).join(sm_bg, suffix = '_bg')
    
    dna_first = dna_first.drop(like="bg$").as_df()
 
    # restore initial order: rna parts comes first
    return dna_first.iloc[:, np.r_[6:11, 0:6, 11]].rename(columns={"Chromosome": "dna_chr",
                                                                   "Start": "dna_bgn",
                                                                   "End": "dna_end"})

def calculate_N2_normalization(rdc_backannot):
    
    """
    Calculates N2 metrics and renormalizes weights to preserve the total number of RNA-DNA contacts

    Parameters
    ----------
    rdc_backannot: pd.DataFrame
        Contacts file with bakground values assigned
        
    Returns
    -------
    pd.DataFrame
    
    """
    
    rdc_backannot["N2_raw"] = 1 / (rdc_backannot["bg_sm"] + 0.5)
    
    library_size = rdc_backannot.shape[0]
    sum_of_weights = sum(rdc_backannot.N2_raw)
    
    rdc_backannot["N2"] = round(rdc_backannot["N2_raw"] * (library_size / sum_of_weights), 3)
    
    rdc_backannot.drop(columns = ["N2_raw"]).to_csv(os.path.join(args.outdir, Path(args.path_rdc).stem + '_N2.tab'), sep = '\t', 
                                                                                                            index=False)
    return

run_normalization(rdc_path = args.path_rdc, 
                  stereogene_path = args.path_Stereogene,
                  config_path = args.path_config,
                  chrsizes_path = args.path_chrsizes,
                  top_n_drop = args.filter_top_n,
                  tail_n_drop = args.filter_tail_n,
                  outdir = args.outdir, stats_outdir = args.stat_folder)
