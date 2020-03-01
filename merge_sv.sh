#!/bin/bash

usage() {
  echo "usage: $0 <in_bed2> <out_bed2>"
  exit 1
}

if [ -z $1 ] || [ -z $2 ];
then
  usage
  exit 1
fi

in_bed2=$1
out_bed2=$2

set -e

raw_bed2=${in_bed2}
raw_bed=${in_bed2}.bed
raw_vcf=${in_bed2}.vcf
merge_vcf=${in_bed2}.merge.vcf
merge_bed=${in_bed2}.merge.bed
merge_file=${in_bed2}.merge.txt
merge_bed2=${out_bed2}

survivor=SURVIVOR
max_dist_between_bp=100
min_supp=2
min_sv_size=1000

cat ${raw_bed2} | awk '{printf("%s\t%s\t%s\n",$1,$2,$4)}' > ${raw_bed}

${survivor} bedtovcf ${raw_bed} DEL ${raw_vcf}
printf "${raw_vcf}\n${raw_vcf}\n" > ${merge_file}
${survivor} merge ${merge_file} ${max_dist_between_bp} ${min_supp} 1 1 0 ${min_sv_size} ${merge_vcf}
python vcf2Bed2.py -v ${merge_vcf} > ${merge_bed2}
