#!/bin/bash



item1=$(echo -e "$1" | cut -f 1)
#root="/mnt/scratcha/fmlab/schnei01/Data/download/syntheticDuplicates/bams/UID-FML-PEO1-FUCCI/"
root="/mnt/scratcha/fmlab/schnei01/Data/download/syntheticDuplicates/bams/UID-FML-PEO1-PLOIDY-2N/"

echo $item1
file1="$root$item1"
#newname="downsample/UID-FML-PEO1-FUCCI-DOWNSAMPLE-50%_$(echo $item1 | cut -f 2,3,4 -d "_")"
newname="downsample/UID-FML-PEO1-PLOIDY-2N-DOWNSAMPLE-25%_$(echo $item1 | cut -f 2,3,4 -d "_")"
echo $file1
echo $newname
# subsample seed is 11
#example: samtools view -s 11.25 $file1 > ${file1}.tmp
#
# for PLOIDY-2N downsample to total reads: 3 754 797
ratio=".25"
samtools view -h -b -o ${newname} -s 11${ratio} $file1
reads=$(samtools view -c $file1)
#ratio=$(echo "scale=2; 3754800/$reads" | bc -l)
ratio=9.0
echo "A"
echo "$(echo "$ratio < 1.00" | bc -l)"
if [ $(bc <<< "$ratio < 1.0") -eq 1 ]; then
	echo "downsampling"
	echo $reads
	echo $ratio
	samtools view -h -b -o ${newname} -s 11${ratio} $file1
else
	# downsample not necessary
	echo "not downsampling"
	echo $reads
	echo $ratio
fi

