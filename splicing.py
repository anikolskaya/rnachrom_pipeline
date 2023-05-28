import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np

input_schema = ['rna_chr',
                'rna_bgn',
                'rna_end',
                'id',
                'rna_strand',
                'rna_cigar',
                'dna_chr',
                'dna_start',
                'dna_end',
                'dna_strand',
                'dna_cigar']

input_dtypes = {'rna_chr': 'category',
                'rna_bgn': 'int',
                'rna_end': 'int',
                'id': 'int',
                'rna_strand': 'category',
                'rna_cigar': 'category',
                'dna_chr': 'category',
                'dna_start': 'int',
                'dna_end': 'int',
                'dna_strand': 'category',
                'dna_cigar': 'category'}


def load_rdc(filename: str) -> pd.DataFrame:
    return pd.read_csv(filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=input_schema,
                       dtype=input_dtypes)

def splicing_cases(df: pd.DataFrame) -> pd.DataFrame:
    df['N_cnt'] = df.rna_cigar.str.count('N')
    df['ID_cnt'] = df.rna_cigar.str.count(r'[ID]')
    return df


def process_simple_splicing(df: pd.DataFrame) -> pd.DataFrame:
    df_ss = df[(df.N_cnt == 1) & (df.ID_cnt == 0)]
    print(df_ss.head())
    df_ss[['N1','N2']] = df_ss['rna_cigar'].str.split(r'[NM]', expand=True).iloc[:,[0,2]]
    df_ss = df_ss.astype({'N1':'int','N2':'int'})
    
    df_ss['rna_end'] = np.select([df_ss.N1 >= df_ss.N2], [df_ss.rna_bgn + df_ss.N1], default = df_ss.rna_end)
    df_ss['rna_bgn'] = np.select([df_ss.N1 < df_ss.N2], [df_ss.rna_end - df_ss.N2], default = df_ss.rna_bgn)
    
    return df_ss


#def stats(df_cnt: pd.DataFrame) -> pd.DataFrame: 
# статистика по сплайсингу

test = load_rdc(sys.argv[1])
test = splicing_cases(test)

match = test[(test.N_cnt == 0) & (test.ID_cnt == 0)]
test_ss = process_simple_splicing(test)
out = pd.concat([match, test_ss]).drop(['N1', 'N2', 'N_cnt', 'ID_cnt'], axis=1)
out.to_csv(Path(sys.argv[1]).stem + '_CIGAR_processed.tsv', sep = '\t', header=False, index=False)
