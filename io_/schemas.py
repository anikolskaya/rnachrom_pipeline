"""
Field names for various genomic tabular files

"""


rdc_BED = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'rna_cigar',
    'dna_chr',
    'dna_start',
    'dna_end',
    'dna_strand',
    'dna_cigar'
]

rdc_dtypes = {
    'rna_chr': 'category',
    'rna_bgn': 'int',
    'rna_end': 'int',
    'id': 'str',
    'rna_strand': 'category',
    'rna_cigar': 'category',
    'dna_chr': 'category',
    'dna_start': 'int',
    'dna_end': 'int',
    'dna_strand': 'category',
    'dna_cigar': 'category'
}

rdc_pyranges_BED = [
    'Chromosome',
    'Start',
    'End',
    'id',
    'Strand',
    'rna_cigar',
    'dna_chr',
    'dna_start',
    'dna_end',
    'dna_strand',
    'dna_cigar']


annot_pyranges_BED = [
    'Chromosome',
    'Start',
    'End',
    'Strand',
    'gene_name',
    'gene_type',
    'source',
]


voted_BED = [
    'rna_chr',
    'rna_bgn',
    'rna_end',
    'id',
    'rna_strand',
    'dna_chr',
    'dna_start',
    'dna_end',
    'dna_strand',
    'gene_name',
    'gene_type' 
]

# UCSC File Format
# https://genome.ucsc.edu/FAQ/FAQformat.html
BED6 = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand"]
