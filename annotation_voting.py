def run_annotation_and_voting(gene_annot_path: str, 
                              contacts_path: str,
                              keep_strand = True,
                              ncpus = 1):
    gene_annot = pr.read_gtf(gene_annot_path)
    cnts = load_rdc(contacts_path)
    
    if keep_strand:
        cnts_annot = annotate_rdc(cnts, gene_annot, cpus = ncpus)
    else:
        cnts_annot = annotate_rdc(cnts, gene_annot, 
                                  strand = 'opposite', cpus = ncpus)
    
    cnts_annot = cnts_annot.as_df()
    no_annot = cnts_annot[cnts_annot['gene_name'] == '-1']
    no_annot.drop(['gene_name', 'gene_length'], axis = 1, inplace = True)
    
    cnts_annot_voted = vote(cnts_annot[cnts_annot['gene_name'] != '-1'])
    
    return cnts_annot_voted, no_annot

def load_rdc(filename: str, 
             BED_columns = ['Chromosome', 'Start', 'End', 'ID', 'Strand']) -> pr.PyRanges:
    rdc = pd.read_csv(filename,
                       header=None,
                       index_col=None,
                       sep='\t',
                       names=input_schema,
                       dtype=input_dtypes).iloc[:,:5]
    rdc.columns = BED_columns
    return pr.PyRanges(rdc)
    

def annotate_rdc(contacts: pr.PyRanges, annot: pr.PyRanges, strand = 'same',
                cpus = 1) -> pr.PyRanges:
    return contacts.join(gene_annot, how = 'left', strandedness = strand,
                         nb_cpu = cpus)[['Strand', 'ID', 'gene_name', 'gene_length']]

def vote(an_contacts: pd.DataFrame) -> pd.DataFrame:
    an_contacts['count'] = an_contacts['gene_name'].map(an_contacts['gene_name'].value_counts())
    an_contacts['density'] = (an_contacts['count'] / 
                                 an_contacts['gene_length']) * 1000000
    
    return rdc_annot_df.sort_values(['ID','density'],ascending=False).groupby('ID').head(1)
