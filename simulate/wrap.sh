#!/bin/sh
int=$1
ref=$2
fastq_1=$3
fastq_2=$4
out=$5
cmd="/y/home/luozhihui/program/bwa-0.7.12/bwa mem -t $int $ref $fastq_1 $fastq_2 >$out"
eval $cmd
