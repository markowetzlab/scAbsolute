#!/bin/bash

DATASET=$1
RESULTPATH=$2
PLOIDY=$3
# DATASET="/home/schnei01/Iapetus/benchmarking/simulations/PEO1-simulated-calls_5_FALSE.csv"
# RESULTPATH="/home/schnei01/Iapetus/benchmarking/results/SCICoNE/"
# PLOIDY=2
echo "Parameters"
echo $DATASET
echo $RESULTPATH
echo $PLOIDY
post=${DATASET##*/}
post=${post%.*}

columns=$(/usr/bin/gawk -F "," '{print NF}' $DATASET | sort -nu | tail -n 1)
rows=$(cat $DATASET | wc -l)
echo $columns
echo $rows

cd $RESULTPATH
echo "Writing output to"
echo "$(pwd -P)"
echo "$post"

/opt/SCICoNE/build/breakpoint_detection --d_matrix_file ${DATASET} --n_bins=$columns --n_cells=$rows --window_size 30 --threshold 3.0 --bp_limit 5000 --bp_min=10 --min_cells=3 --verbosity=1 --evaluate_peaks=True --postfix="$post"

echo "${post}_segmented_region_sizes.txt"
python /opt/SCICoNE/scripts/segment_counts.py ${DATASET} "$RESULTPATH/${post}_segmented_region_sizes.txt"
mv -f ${DATASET%.csv}_segmented_counts.txt $RESULTPATH/

echo "DEBUG inference"
echo "${post}_segmented_regions.txt"
echo $(cat "${post}_segmented_region_sizes.txt" | wc -l)
echo ${DATASET}
echo $PLOIDY
echo "${RESULTPATH}/${post}_segmented_counts.txt"
echo "$(cat "${post}_segmented_regions.txt" | wc -l)"
echo "START"

/opt/SCICoNE/build/inference --n_cells $rows  --n_regions $(cat "${post}_segmented_region_sizes.txt" | wc -l) --n_iters 10000 --n_nodes 100 --ploidy $PLOIDY --verbosity 2 --seed 42 --copy_number_limit 8 --d_matrix_file $(echo "${RESULTPATH}/${post}_segmented_counts.txt" | tr -s / ) --region_sizes_file="${post}_segmented_region_sizes.txt" --postfix="$post"

mv -f "${post}_inferred_cnvs.csv" "${post}.scicone.csv"
mv -f "${post}_cell_node_ids.tsv" "${post}.scicone.nodes.csv"
mv -f "${post}_tree_inferred.txt" "${post}.scicone.tree"

echo "COMPLETED"
