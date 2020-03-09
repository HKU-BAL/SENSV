#!/bin/bash

usage() {
  echo "usage: $0 <read_name> <fastq>"
  exit 1
}

if [ -z $1 ] || [ -z $2 ];
then
  usage
  exit 1
fi

readname=$1
fastq=$2

num=`look -b ${readname} ${fastq}.idx|head -1l|cut -d' ' -f2`
grabix grab ${fastq} $num $((num+3))
