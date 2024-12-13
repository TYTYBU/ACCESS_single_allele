#!/bin/bash -l

#$ -N cram2fastq
#$ -j y
#$ -V
#$ -q all.q@compute-0-0,all.q@compute-0-1
#$ -pe serial 64
#$ -M tyu7@bwh.harvard.edu
#$ -wd /net/bgm/sherwood/
#$ -o /net/bgm/sherwood/job_out

DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
mkdir -p $DATA_PATH/fastq/fastqc


# for f in $DATA_PATH/*.ucram
# do
# 	f2=$(basename $f .ucram)
# 	samtools view --tag a3 $f | samtools fastq | gzip -f > $DATA_PATH/fastq/${f2}.full_len.fastq.gz &
# 	# fastqc -o $DATA_PATH/fastq/fastqc/ $DATA_PATH/fastq/${f2}.fastq.gz
# done
# wait

# for f in $DATA_PATH/*.ucram
# do
# 	f2=$(basename $f .ucram)
# 	samtools fastq $f | gzip -f > $DATA_PATH/fastq/${f2}.fastq.gz &
# done
# wait


# SAMPLE_SHEET=$DATA_PATH/N5_index_to_names.txt
# mkdir -p $DATA_PATH/fastq/all
# mkdir -p $DATA_PATH/fastq/full_len

# [ ! -f $SAMPLE_SHEET ] && { echo "$SAMPLE_SHEET file not found"; exit 99; }
# while IFS=',' read -r sample_name index_seq
# do
# 	echo "sample name: $sample_name"
# 	cat $DATA_PATH/fastq/*${index_seq}.trimmed.fastq.gz > $DATA_PATH/fastq/all/${sample_name}.fastq.gz &
# 	cat $DATA_PATH/fastq/*${index_seq}.trimmed.full_len.fastq.gz > $DATA_PATH/fastq/full_len/${sample_name}.full_len.fastq.gz &
# done < <(cat $SAMPLE_SHEET)
# wait

for f in $DATA_PATH/fastq/all/*.fastq.gz
do
	echo $f $(( $(zcat $f | wc -l) / 4 )) &
done
wait

# for f in $DATA_PATH/fastq/full_len/*.full_len.fastq.gz
# do
# 	echo $f $(( $(zcat $f | wc -l) / 4 )) &
# done
# wait



