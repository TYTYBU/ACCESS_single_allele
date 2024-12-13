#!/bin/bash -l

#$ -N dedup
#$ -j y
#$ -V
#$ -q all.q@compute-0-0,all.q@compute-0-1
#$ -pe serial 64
#$ -M tyu7@bwh.harvard.edu
#$ -wd /net/bgm/sherwood/
#$ -o /net/bgm/sherwood/job_out

DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
mkdir -p $DATA_PATH/dedup_fastq/all
mkdir -p $DATA_PATH/dedup_fastq/full_len
mkdir -p $DATA_PATH/temp


# for f in $DATA_PATH/fastq/all/*.fastq.gz
# do
# 	f2=$(basename $f .fastq.gz)
# 	echo $f2
# 	zcat $f > $DATA_PATH/temp/temp.fastq
# 	czid-dedup --inputs $DATA_PATH/temp/temp.fastq -o $DATA_PATH/dedup_fastq/all/${f2}.dedup.fastq
# 	# gzip -f $DATA_PATH/dedup_fastq/all/${f2}.dedup.fastq
# 	pigz -f -p 64 $DATA_PATH/dedup_fastq/full_len/${f2}.dedup.fastq
# done

for f in $DATA_PATH/fastq/full_len/*.fastq.gz
do
	f2=$(basename $f .fastq.gz)
	echo $f2
	zcat $f > $DATA_PATH/temp/temp.fastq
	czid-dedup --inputs $DATA_PATH/temp/temp.fastq -o $DATA_PATH/dedup_fastq/full_len/${f2}.dedup.fastq
	pigz -f -p 64 $DATA_PATH/dedup_fastq/full_len/${f2}.dedup.fastq
done




