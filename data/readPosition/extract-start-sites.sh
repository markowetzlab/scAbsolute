#!/bin/bash

#A fairly simple script which takes an indexed, co-ordinate sorted bam file as input, marks and remove duplicates, and then produces a table of reads including the reference and genomic start site for each read

# Original author: Thomas Bradley
# Modified by: Michael Schneider 

input_bam=$1
relevant_reference_sequences=$2
workdir=$(dirname "$0")
#echo $workdir

# mark duplicates
picard MarkDuplicates I=${input_bam} O=${input_bam}.deduplicated M=${input_bam}.duplicate_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
# java -jar ./picard.jar MarkDuplicates I=${input_bam} O=${input_bam}.deduplicated M=${input_bam}.duplicate_metrics.txt REMOVE_DUPLICATES=true

# convert chromosome/contig co-ordinates to genome co-ordinates
samtools view -H ${input_bam}.deduplicated | grep '^@SQ' | awk '{FS="\t"}{OFS="\t"}{print $2,$3}'  > ${input_bam}.reference_sizes.txt

#
if [ $(samtools view -c -f 1 ${input_bam}) = "0" ]; then
    # single-end reads
    echo "single-end reads"
    # get table from deduplicated bam
    # we use read length (length($10)) as opposed to fragment size in paired end sequencing
    samtools view -F "0x400" ${input_bam}.deduplicated | awk '{FS="\t"}{OFS="\t"}{print $1,$3,$4,length($10)}' > ${input_bam}.deduplicated.read_start_sites.txt

else
    # paired-end reads
    echo "paired-end reads"
    # get table from deduplicated bam
    # samtools view ${input_bam}.deduplicated
    samtools view -F "0x400" -f 67 ${input_bam}.deduplicated | awk '{FS="\t"}{OFS="\t"}{print $1,$3,$4,$9}' > ${input_bam}.deduplicated.read_start_sites.txt
fi

# convert read start positions into genomic co-ordinates
Rscript --vanilla $workdir/get_table.R ${input_bam}.reference_sizes.txt ${relevant_reference_sequences} ${input_bam}.deduplicated.read_start_sites.txt ${input_bam%.bam}.position.tsv

# clean up
rm -f ${input_bam}.deduplicated
rm -f ${input_bam}.duplicate_metrics.txt
rm -f ${input_bam}.deduplicated.read_start_sites.txt
rm -f ${input_bam}.reference_sizes.txt

echo "COMPLETED"
