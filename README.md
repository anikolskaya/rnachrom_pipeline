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
