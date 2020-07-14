#!/bin/sh
coverage=$1
purity=$2
EAGLE=/y/home/fanxp/software/EAGLE
date
echo "Simulation for coverage: ${coverage}, purity: ${purity}"
cmd="mkdir purity_$purity"
echo "run $cmd"
eval $cmd
cmd="cd purity_$purity"
echo "run $cmd"
eval $cmd

# calcualte tumor coverage and normal coverage
tumor_coverage=$(echo $coverage $purity | awk '{ printf "%0.2f\n", $1*$2}')
normal_coverage=$(echo $coverage $purity | awk '{ printf "%0.2f\n", (1-$2)*$1}')
echo "tumor sample coverage: ${tumor_coverage}. normal sample coverage: ${normal_coverage}"

#create refence genome and normal cell
cmd="${EAGLE}/bin/configureEAGLE.pl --run-info=${EAGLE}/share/EAGLE/RunInfo/RunInfo_PairedReadsBarcode8x32Tiles.xml --reference-genome=../hs37d5_namechr.fa --variant-list=../total_vcf.vcf --coverage-depth=${normal_coverage} --genome-mutator-options=\"--organism-ploidy=2 --ploidy-chromosome=chrY --ploidy-level=0\" EAGLE_normal"
echo "run $cmd"
eval $cmd
cmd="cd EAGLE_normal"
echo "run $cmd"
eval $cmd
cmd="make sample"
echo "run $cmd"
eval $cmd
cmd="make bam -j 30"
echo "run $cmd"
eval $cmd

# generate CNA.vcf
cmd="cd .."
echo "run $cmd"
eval $cmd
cmd="${EAGLE}/libexec/EAGLE/applyCopyNumber.pl -i ../singleclone_CN.tab -o CNA.vcf"
echo "run $cmd"
eval $cmd

#create tumor cell
cmd="${EAGLE}/bin/configureEAGLE.pl --run-info=${EAGLE}/share/EAGLE/RunInfo/RunInfo_PairedReadsBarcode8x32Tiles.xml --reference-genome=./EAGLE_normal/sample_genome --variant-list=./CNA.vcf --coverage-depth=${tumor_coverage} --genome-mutator-options=\"--organism-ploidy=1 --ploidy-chromosome=chrY --ploidy-level=0\" EAGLE_tumor"
echo "run $cmd"
eval $cmd
cmd="cd EAGLE_tumor"
echo "run $cmd"
eval $cmd
cmd="make bam -j 30"
echo "run $cmd"
eval $cmd

cmd="cd .."
echo "run $cmd"
eval $cmd

#use picard to split bam to p1 and p2
cmd="java -jar /y/home/fanxp/software/picard/picard.jar SamToFastq I=EAGLE_normal/eagle.bam FASTQ=normal_eagle_1.fastq SECOND_END_FASTQ=normal_eagle_2.fastq VALIDATION_STRINGENCY=LENIENT"
echo "run $cmd"
eval $cmd
cmd="java -jar /y/home/fanxp/software/picard/picard.jar SamToFastq I=EAGLE_tumor/eagle.bam FASTQ=tumor_eagle_1.fastq SECOND_END_FASTQ=tumor_eagle_2.fastq VALIDATION_STRINGENCY=LENIENT"
echo "run $cmd"
eval $cmd

#add prefix to distinct tumor cell and normal cell 
cmd="/y/Sunset/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/addPrefix.pl normal_eagle_1.fastq normal_eagle_1_add.fastq 1"
echo "run $cmd"
eval $cmd
cmd="/y/Sunset/tcga/gdc-tools/Puricise/simulation_pipeline/sbin/addPrefix.pl normal_eagle_2.fastq normal_eagle_2_add.fastq 1"
echo "run $cmd"
eval $cmd

#mixture tumor cell and normal cell to generate tumor sample
cmd="cat normal_eagle_1_add.fastq tumor_eagle_1.fastq > mix_1.fastq"
echo "run $cmd"
eval $cmd
cmd="cat normal_eagle_2_add.fastq tumor_eagle_2.fastq > mix_2.fastq"
echo "run $cmd"
eval $cmd

#run bwa
cmd="bwa mem -t 10 /y/Sunset/db/individual_sequence/1_hs_genome_hs37d5.fa mix_1.fastq mix_2.fastq > tum_aligned.sam"
echo "run $cmd"
eval $cmd

#convert sam to bam,sort and index
cmd="samtools view -bS tum_aligned.sam -o tum_aligned.bam"
echo "run $cmd"
eval $cmd
cmd="samtools sort tum_aligned.bam -@ 20 -m 2G -o tum_aligned.sort.bam"
echo "run $cmd"
eval $cmd
cmd="samtools index tum_aligned.sort.bam"
echo "run $cmd"
eval $cmd
date
