Excluded regions:
https://github.com/hall-lab/speedseq
https://github.com/hall-lab/speedseq/blob/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed
exclude regions, consider only regions > 3kb

`awk -v OFS="\t" '{$4=$3-$2}1' ceph18.b37.lumpy.exclude.2014-01-15.bed.txt | sort -nr -k 4 | head -n 1120 > excluded_regions.bed`


Gencode annotations:
source:
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz

cat gencode.v36.annotation.gff3 | grep -P "gene\t" | grep -v "^chrM" | cut -f 1,2,4,5,9 | tr ";" "\t" | cut -f 1,2,3,4,6,8 | sed -e 's/gene_id=//g' -e 's/gene_name=//g' >! gencode.v36.tsv
cat gencode.v19.annotation.gff3 | grep -P "gene\t" | grep -v "^chrM" | cut -f 1,2,4,5,9 | tr ";" "\t" | cut -f 1,2,3,4,6,10 | sed -e 's/gene_id=//g' -e 's/gene_name=//g' >! gencode.v19.tsv
