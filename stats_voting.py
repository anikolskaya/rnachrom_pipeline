import pyranges as pr
import pandas as pd
import pathlib
import os
from natsort import natsort_keygen

import pathlib
import os
from natsort import natsort_keygen

def calculate_stats(frames_d: dict,
                    outdir_path = './stats_2'):
    """
    Calculate statistics for annotation and voting processes.
    

    Parameters
    ----------
    frames_d: dict
        Dict from run_annotation_and_voting command
        
    outdir_path:
        Path to the output directory
    
    """
    
    pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)
    
    frames_d['voted_annot'].sort_values('gene_name', ascending=True, inplace = True)
    
    stats_by_gene_type(frames_d['voted_annot']).to_csv(os.path.join(outdir_path,'counts_gene_type.tab'),
                                                       index = True, 
                                                       sep = '\t')
    
    density_by_gene(frames_d['voted_annot']).to_csv(os.path.join(outdir_path,'genes_density.tab'), 
                                                    index = False, 
                                                    sep = '\t')
    
    counts_by_gene(frames_d['voted_annot']).to_csv(os.path.join(outdir_path,'counts_by_genes.tab'), 
                                                   index = True, 
                                                   sep = '\t')
    
    counts_by_chr(frames_d).to_csv(os.path.join(outdir_path,'counts_by_chr.tab'), 
                                               index = True, 
                                               sep = '\t')
    return
    



def stats_by_gene_type(vote_sorted: pd.DataFrame):
    
    return vote_sorted.groupby(['gene_type']).agg(n_counts=('gene_name', 'count'), 
                                           n_genes=('gene_name', 'nunique'))
                                                                                                                                                            
    


def density_by_gene(vote_sorted: pd.DataFrame):
    
    return vote_sorted[['Chromosome', 'gene_name', 'gene_type', 'density']].drop_duplicates('gene_name')   
    


def counts_by_gene(vote_sorted: pd.DataFrame, filename = 'counts_by_genes.tab'):
    
    return vote_sorted.groupby('gene_name').agg(count=('gene_name', 'count'))     


def counts_by_chr(d: dict):
    d.pop('voted_annot', None)
    all_frames_stat = dict(map(lambda item: (item[0], item[1].groupby('Chromosome')['id'].nunique()), d.items()))
    
    return pd.DataFrame.from_dict(all_frames_stat).rename(columns={'complement_annot': 'Nannot_wrong_strand',
                                                                   'selected_annot': 'Nannot_correct_strand',
                                                                   'no_annot': 'Nnoannot',
                                                                   'voted_noribo': 'NnoRibo'}).sort_values(by="Chromosome", key=natsort_keygen())
    