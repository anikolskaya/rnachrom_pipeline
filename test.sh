
#python3 splicing.py -h
#python3 blacklisting.py -h
#python3 annotation_voting.py -h
#python3 splicing.py -h

python3 splicing.py test_files/SRR9201799.tab
python3 blacklisting.py ./results/SRR9201799_CIGAR_processed.tsv hg38
python3 -W ignore annotation_voting.py ./results/SRR9201799_CIGAR_processed_blacklst.tsv test_files/gencode.v43.annotation.gene.gtf --annot_file_format GTF


python3 normalize_N2.py ./results/GRID_mES_pipe.chr19_voted_annot.tab /mnt/lustre/suvorova/tools/src/Smoother ./test_files/run_stereogene/cfg.cfg ./test_files/run_stereogene/hg38_chrsize.tab --filter_top_n 1 --filter_tail_n 1
python3 scaling.py rdc_norm_path chrsizes --threads 4
