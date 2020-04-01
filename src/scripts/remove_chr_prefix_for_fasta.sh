#!/bin/bash

# usage example ./remove_chr_prefix_for_fasta.sh "<fasta file absolute path>" ./GRCh38
usage() {
  echo "usage: $0 <fasta_file_path> <out_prefix>"
  exit 1
}

if [ -z $1 ] || [ -z $2 ];
then
  usage
  exit 1
fi

set -e

fasta_file_path=$1
out_prefix=$2

cat ${fasta_file_path} | sed 's/>chr/>/g' > ${out_prefix}.no_chr.fa
samtools faidx ${out_prefix}.no_chr.fa
