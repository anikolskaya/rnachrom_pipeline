
#python3 splicing.py -h
#python3 blacklisting.py -h
#python3 annotation_voting.py -h

python3 splicing.py test_files/SRR9201799.tab
python3 blacklisting.py ./results/SRR9201799_CIGAR_processed.tsv hg38
python3 -W ignore annotation_voting.py ./results/SRR9201799_CIGAR_processed_blacklst.tsv test_files/gencode.v43.annotation.gene.gtf --annot_file_format GTF
