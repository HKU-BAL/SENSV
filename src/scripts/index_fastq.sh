#!/bin/bash

usage() {
  echo "usage: $0 <in_fastq.gz> <out_prefix> <nthread>"
  exit 1
}

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ];
then
  usage
  exit 1
fi

set -e

in_fastq_gz=$1
out_prefix=$2
nthread=$3

out_gz=${out_prefix}.fastq.gz

echo "out_gz=$out_gz"

pigz -fdc -p ${nthread} ${in_fastq_gz} | \
tee >(bgzip -@ ${nthread} -c > ${out_gz}) | \
awk 'NR%4==1 {print substr($1,2),NR}' | \
sort --parallel=${nthread} > ${out_gz}.idx

grabix index ${out_gz}
