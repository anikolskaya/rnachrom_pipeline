# rnachrom_pipeline
usage: splicing.py [-h] [--outdir [OUTDIR]] path_rdc_primary  
Parse CIGAR field and process spliced contacts  
positional arguments:  
  path_rdc_primary   Path to the primary rna-dna contacts file. 
  
  
usage: blacklisting.py [-h] [--outdir [OUTDIR]] path_rdc genome

Removes contacts in which DNA parts overlap with ENCODE blacklisted regions

positional arguments:
  path_rdc           Path to the rna-dna contacts file
  genome             Genome abbreviation (ENCODE blacklisted regions are
                     available for hg19, hg38, mm9, mm10)

optional arguments:
  -h, --help         show this help message and exit
  --outdir [OUTDIR]  A folder in which to store the results
  
  
usage: annotation_voting.py [-h] [--annot_file_format [ANNOT_FILE_FORMAT]] [--outdir [OUTDIR]] rdc gene_annot  
Annotate RNA-parts of contacts, perform a voting procedure and deplete contacts   
positional arguments:  
  rdc                   Path to rna-dna contacts file  
  gene_annot            Path to gene annotation  
optional arguments:  
  -h, --help            show this help message and exit. 
  --annot_file_format [ANNOT_FILE_FORMAT]. 
                        Gene annotation file format: GTF or BED (default: GTP)  
  --outdir [OUTDIR]     A folder in which to store results  

usage: normalize_N2.py [-h] [--filter_top_n [FILTER_TOP_N]] [--filter_tail_n [FILTER_TAIL_N]] [--outdir [OUTDIR]]
                       [--stat_folder [STAT_FOLDER]]
                       path_rdc path_Stereogene path_config path_chrsizes

Calculate N2 background normalization

positional arguments:
  path_rdc              Path to the rna-dna contacts file
  path_Stereogene       Path to the Stereogene smoother script
  path_config           Path to the config template file
  path_chrsizes         Path to the file with chromosome sizes

optional arguments:
  -h, --help            show this help message and exit
  --filter_top_n [FILTER_TOP_N]
                        Drop the most contacting RNAs (rank, default: 50)
  --filter_tail_n [FILTER_TAIL_N]
                        Drop the least contacting RNAs (rank, default: 1000)
  --outdir [OUTDIR]     A folder in which to store the results
  --stat_folder [STAT_FOLDER]
                        A folder in which to store stats

usage: scaling.py [-h] [--factor [FACTOR]] [--threads [THREADS]]
                  [--outdir [OUTDIR]]
                  path_rdc_norm path_chrsizes

Calculates inner- and inter-chromosomal scaling based on mRNA contacts density
decay

positional arguments:
  path_rdc_norm        Path to the rna-dna contacts file,
                       background_normalized
  path_chrsizes        Path to the file with chromosome sizes

optional arguments:
  -h, --help           show this help message and exit
  --factor [FACTOR]    Scale factor
  --threads [THREADS]
  --outdir [OUTDIR]    A folder in which to save results
