#!/bin/bash

echo $1

#. /home/schnei01/.miniconda3/etc/profile.d/conda.sh
#eval "$(/home/schnei01/.miniconda3/bin/conda shell.bash hook)"
#/home/schnei01/.miniconda3/bin/conda activate chisel
spack unload py-matplotlib@3.5.1
spack unload py-numpy@1.22.0
spack unload py-cycler@0.11.0
spack unload py-kiwisolver@1.3.2
spack unload py-pillow@8.4.0
spack unload py-packaging@21.3
spack unload py-fonttools@4.28.1
spack unload py-pyparsing@3.0.6
spack unload py-python-dateutil@2.8.2
spack unload py-setuptools@57.4.0
spack unload py-six@1.16.0
spack unload python@3.9.9

#/usr/bin/awk '{gsub(/^chr/,""); print}' LIFTOVER/$1.mpileup.phased.liftover.vcf > CHISEL/$1/${1%.vcf}.nochr.vcf
#rsync -aP HRC/$1.mpileup.HRC-Sanger.vcf.gz CHISEL2/$1/

mkdir -p CHISEL/$1
cd CHISEL/$1

chisel_nonormal -t /mnt/scratchb/fmlab/schnei01/Data/project1/CHISEL/ACT/NONORMAL2/bams/$1.bam -r /mnt/scratchb/fmlab/schnei01/Data/reference/GRCh37/GRCh37.fa --seed 2022 -l /mnt/scratchb/fmlab/schnei01/Data/project1/CHISEL/ACT/NONORMAL2/TOPMED/$1.merged.phased.liftover.nochr.tsv --jobs 11 --chromosomes "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
