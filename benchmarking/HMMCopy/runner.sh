 single_cell hmmcopy \
 --input_yaml UID-NNA-TN5.yaml \
 --library_id NAVINACT --maxjobs 24 --nocleanup \
 --sentinel_only --loglevel DEBUG  --submit local \
 --tmpdir tmp --pipelinedir pipeline --out_dir output/UID-NNA-TN5/ \
 --config_override '{"refdir": "/mnt/scratchb/fmlab/schnei01/Data/project1/HMMCopy/refdata/",  "hmmcopy": {"gc_wig_file": "/mnt/scratchb/fmlab/schnei01/Data/project1/HMMCopy/refdata/human/GRCh37.gc.ws_500000.wig", "map_wig_file": "/mnt/scratchb/fmlab/schnei01/Data/project1/HMMCopy/refdata/human/GRCh37.map.ws_125_to_500000.wig", "ref_genome": "/mnt/scratchb/fmlab/schnei01/Data/project1/HMMCopy/refdata/human/GRCh37.fa"}}'
#, "chromosomes": ["1"]}}'
# --tmpdir tmp --pipelinedir pipeline --output_prefix output/UID-DLP-SA922/ \
#
