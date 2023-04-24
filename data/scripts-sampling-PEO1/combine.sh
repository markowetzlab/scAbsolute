#!/bin/bash

while read -r line; do
	#echo $line;

	item1=$(echo -e "$line" | cut -f 1)
	item2=$(echo -e "$line" | cut -f 2)
	item3=$(echo -e "$line" | cut -f 3)
	root="/mnt/scratcha/fmlab/schnei01/Data/download/syntheticDuplicates/bams/UID-FML-PEO1-FUCCI/"

	echo $item1
	echo $item2
	file1="$root$item1"
	file2="$root$item2"
	newname="synthetic/UID-FML-PEO1-SYN-2N-0.50_SLX-00000_${item3}_$(echo $item1 | cut -f 2 -d "/" | tr _ \* | sed 's/\.[^.]*$//')-*-$(echo $item2 | cut -f 2 -d "/" | tr _ \* | sed 's/\.[^.]*$//').bam"
	#echo $file1
	#echo $newname
	# subsample seed is 11
	#example: samtools view -s 11.25 $file1 > ${file1}.tmp

	samtools view -s 11.50 -bo "${file1}.tmp" $file1;
	samtools view -s 11.50 -bo "${file2}.tmp" $file2;
	sleep 3
	samtools merge -o $newname ${file1}.tmp ${file2}.tmp
	sleep 3;
	rm -f ${file1}.tmp
	rm -f ${file2}.tmp

	echo ""
done < $1
