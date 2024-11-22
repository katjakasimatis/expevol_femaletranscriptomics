#!/bin/bash

##REMOVE ADAPTERS WITH SKEWER##
module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 skewer/0.2.2

dir=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics/5909
out_dir=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics/5909
final_dir=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics/trimmedReads2

#check adapters are in files
grep -c "AGATCGGAAGAG" $dir/*.fastq

#make an array for each file in the directory "dir" that ends in _R1.fastq
files=(${dir}/*.fastq)

#pulls out one element at a time of the array
file=${files[${SLURM_ARRAY_TASK_ID}]}

name=`basename ${file} *_R1_001.fastq`

#remove TrueSeq adapters + qc
skewer -x AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -y AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -t 12 -l 30 -r 0.01 -d 0.01 -q 10 -o ${out_dir}/${name} $file

#round 2
files2=(${out_dir}/*-trimmed.fastq)
file2=${files2[${SLURM_ARRAY_TASK_ID}]}
name2=`basename ${file2} *-trimmed.fastq`

skewer -x GATCGGAAGAGCACACGTCTGAACTCCAGTCACAACAGGTGATCTCGTAT -t 12 -l 30 -r 0.01 -d 0.01 -q 10 -o ${final_dir}/${name2} $file2




##FAST-QC##
module load java fastqc

LISTFILES=(trimmedReads2/*.fastq) 

file=${LISTFILES[$SLURM_ARRAY_TASK_ID]}

fastqc -o fastqc_out/trimmed_data $file




##MAP READS USING STAR##
module load easybuild icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132 STAR/2.5.3a

dirG=/projects/phillipslab/shared/kasimatis_expevol_genomics
dirT=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics

#Create a genome index using STAR
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $dirT/Ce_star_ref1/ --genomeFastaFiles $dirG/ref_genom_WS274/c_elegans.PRJNA13758.WS274.genomic.fa

#make an array for each file in the directory "trimmedReads" that ends in -trimmed-pair1.fastq
files=(${dirT}/trimmedReads2/*-trimmed-pair1.fastq)

#pulls out one element at a time of the array
file=${files[${SLURM_ARRAY_TASK_ID}]}

name=`basename ${file} -trimmed-pair1.fastq`

#Map reads: pass 1
cd $dirT/star_pass1
STAR --runThreadN 40 --genomeDir $dirT/Ce_star_ref1/ --readFilesIn ${file} ${file/trimmed-pair1.fastq/trimmed-pair2.fastq} --outFileNamePrefix ${name} --outSAMtype BAM SortedByCoordinate

#Make a new genome
cd $dirT/star_pass2
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $dirT/Ce_star_ref2/ --genomeFastaFiles $dirG/ref_genom_WS274/c_elegans.PRJNA13758.WS274.genomic.fa --sjdbFileChrStartEnd $dirT/star_pass1/*.tab

#Map reads: pass 2
cd $dirT/star_pass2
STAR --runThreadN 40 --genomeDir $dirT/Ce_star_ref2 --readFilesIn ${file} ${file/trimmed-pair1.fastq/trimmed-pair2.fastq} --outFileNamePrefix ${name} --outSAMtype BAM SortedByCoordinate




##CREATE TABLE OF GENE COUNTS USING HTSEQ##
module load easybuild icc/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132
module load HTSeq/0.9.1-Python-3.6.1

dir=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics/star_pass2
dirG=/projects/phillipslab/shared/kasimatis_expevol_genomics
out_dir=/projects/phillipslab/shared/kasimatis_expevol_transcriptomics/htseq_raw_counts

#make an array for each file in the directory "trimmedReads" that ends in -trimmed-pair1.fastq
files=(${dir}/*.out.bam)

#pulls out one element at a time of the array
file=${files[${SLURM_ARRAY_TASK_ID}]}

name=`basename ${file} .sortedByCoord.out.bam`

htseq-count -s no -r pos -i ID -t gene -f bam ${dir}/*.out.bam ${dirG}/ref_genom_WS274/c_elegans_LongestIsoform_annotations.gff3 > ${out_dir}/20220617_Gene_table_counts.txt




