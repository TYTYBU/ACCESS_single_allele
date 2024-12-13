#!/bin/bash -l

#$ -N ATACseq_iter
#$ -j y
#$ -V
#$ -q all.q@compute-0-0,all.q@compute-0-1
#$ -pe serial 96
#$ -m ae
#$ -M tyu7@bwh.harvard.edu
#$ -wd /net/bgm/sherwood/NGS_analysis_proj
#$ -o /net/bgm/sherwood/job_out

eval "$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)"
conda activate cutadapt

PROJ_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
SCRIPT_PATH=/net/bgm/sherwood/NGS_analysis_proj/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
GENOME_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/genome
BOWTIE2_INDEX_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/bowtie2_index
NUM_THREADS=96


DATA_PATH=$PROJ_PATH
RAW_PATH=$DATA_PATH/dedup_fastq/all
MERGED_BAM_OUT=$DATA_PATH/bowtie2_out/bam/merged/all


## bowtie2 iterative mapping 
mkdir -p $DATA_PATH/bowtie2_out/bam/mapped_only
mkdir -p $DATA_PATH/bowtie2_out/bam/merged/all
mkdir -p $DATA_PATH/bowtie2_out/fastq
for f in $RAW_PATH/*.fastq.gz # change path depending on whether cutadapt is necessary
do
	f2=$(basename $f .fastq.gz)

	sam_prefix=${f2}.hg38_core
	bam_al=${f2}.hg38_core.al.bam
	out_fq_suffix=hg38_core.un
	bowtie2 -p $(($NUM_THREADS / 3)) --very-sensitive -x $BOWTIE2_INDEX_PATH/hg38_core -U $f | samtools view -@ $(($NUM_THREADS / 3)) -bS - | samtools sort -@ $(($NUM_THREADS / 3)) - -o $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam
	samtools fastq -@ $NUM_THREADS -N -f 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/fastq/${f2}.${out_fq_suffix}.fastq
	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/bam/mapped_only/${bam_al}

	in_fq_suffix=$out_fq_suffix
	sam_prefix=${f2}.un.mp5.hg38_core
	bam_al=${f2}.hg38_core.al2.bam
	out_fq_suffix=hg38_core.un2
	bowtie2 -p $(($NUM_THREADS / 3)) --very-sensitive --mp 5 -x $BOWTIE2_INDEX_PATH/hg38_core -U $DATA_PATH/bowtie2_out/fastq/${f2}.${in_fq_suffix}.fastq | samtools view -@ $(($NUM_THREADS / 3)) -bS - | samtools sort -@ $(($NUM_THREADS / 3)) - -o $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam
	samtools fastq -@ $NUM_THREADS -N -f 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/fastq/${f2}.${out_fq_suffix}.fastq
	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/bam/mapped_only/${bam_al}

	in_fq_suffix=$out_fq_suffix
	sam_prefix=${f2}.un2.mp4.hg38_core
	bam_al=${f2}.hg38_core.al3.bam
	out_fq_suffix=hg38_core.un3
	bowtie2 -p $(($NUM_THREADS / 3)) --very-sensitive --mp 4 -x $BOWTIE2_INDEX_PATH/hg38_core -U $DATA_PATH/bowtie2_out/fastq/${f2}.${in_fq_suffix}.fastq | samtools view -@ $(($NUM_THREADS / 3)) -bS - | samtools sort -@ $(($NUM_THREADS / 3)) - -o $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam
	samtools fastq -@ $NUM_THREADS -N -f 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/fastq/${f2}.${out_fq_suffix}.fastq
	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/bam/mapped_only/${bam_al}

	in_fq_suffix=$out_fq_suffix
	sam_prefix=${f2}.un3.mp3.hg38_core
	bam_al=${f2}.hg38_core.al4.bam
	out_fq_suffix=hg38_core.un4
	bowtie2 -p $(($NUM_THREADS / 3)) --very-sensitive --mp 3 -x $BOWTIE2_INDEX_PATH/hg38_core -U $DATA_PATH/bowtie2_out/fastq/${f2}.${in_fq_suffix}.fastq | samtools view -@ $(($NUM_THREADS / 3)) -bS - | samtools sort -@ $(($NUM_THREADS / 3)) - -o $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam
	samtools fastq -@ $NUM_THREADS -N -f 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/fastq/${f2}.${out_fq_suffix}.fastq
	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/bam/mapped_only/${bam_al}

	in_fq_suffix=$out_fq_suffix
	sam_prefix=${f2}.un4.mp2.hg38_core
	bam_al=${f2}.hg38_core.al5.bam
	out_fq_suffix=hg38_core.un5
	bowtie2 -p $(($NUM_THREADS / 3)) --very-sensitive --mp 2 -x $BOWTIE2_INDEX_PATH/hg38_core -U $DATA_PATH/bowtie2_out/fastq/${f2}.${in_fq_suffix}.fastq | samtools view -@ $(($NUM_THREADS / 3)) -bS - | samtools sort -@ $(($NUM_THREADS / 3)) - -o $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam
	samtools fastq -@ $NUM_THREADS -N -f 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/fastq/${f2}.${out_fq_suffix}.fastq
	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bowtie2_out/bam/${sam_prefix}.bam > $DATA_PATH/bowtie2_out/bam/mapped_only/${bam_al}

	samtools merge -@ $NUM_THREADS -f $MERGED_BAM_OUT/${f2}.hg38_core.merged.bam $DATA_PATH/bowtie2_out/bam/mapped_only/${f2}.hg38_core.al.bam $DATA_PATH/bowtie2_out/bam/mapped_only/${f2}.hg38_core.al2.bam $DATA_PATH/bowtie2_out/bam/mapped_only/${f2}.hg38_core.al3.bam $DATA_PATH/bowtie2_out/bam/mapped_only/${f2}.hg38_core.al4.bam $DATA_PATH/bowtie2_out/bam/mapped_only/${f2}.hg38_core.al5.bam
done

for f in $RAW_PATH/*.fastq.gz
do
	f2=$(basename $f .fastq.gz)
	samtools index $MERGED_BAM_OUT/${f2}.hg38_core.merged.bam &
done
wait

md5sum $MERGED_BAM_OUT/*.bam
rm -f $DATA_PATH/bowtie2_out/bam/mapped_only/*
rm -f $DATA_PATH/bowtie2_out/bam/*.bam
rm -f $DATA_PATH/bowtie2_out/fastq/*








# # for f in $DATA_PATH/bowtie2_out/bam/merged/*.merged.bam
# # do
# # 	f2=$(basename $f .merged.bam)
# # 	samtools index -@ $NUM_THREADS $f
# # 	samtools view -@ $NUM_THREADS -b --regions-file $BED_PATH/combined_TF_flank500nt.Ellie_ACCESS_012524.bed $f > $DATA_PATH/bowtie2_out/bam/merged/${f2}.combined_TF_flank500nt.bam
# # done

# mkdir -p $DATA_PATH/bowtie2_out/baseCt/merged
# for f in $DATA_PATH/bowtie2_out/bam/merged/*.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr_v2.pl $f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa $DATA_PATH/bowtie2_out/baseCt/merged/${f2} &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/baseCt/merged/*.baseCt.txt
# do
# 	gzip -f $f &
# done
# wait











# # if need tp separate WGS and ATAC samples
# mkdir -p $DATA_PATH/bowtie2_out/mpileup_out/ATACseq
# for f in $DATA_PATH/bowtie2_out/bam/merged/*ATAC*.merged_filtered.bam
# do
# 	f2=$(basename $f .hg38_core.merged_filtered.bam)
# 	samtools mpileup --ff UNMAP -B -Q 0 -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/${f2}.mpileup_out.txt &
# done

# mkdir -p $DATA_PATH/bowtie2_out/mpileup_out/WGS
# for f in $DATA_PATH/bowtie2_out/bam/merged/*WGS*.merged_filtered.bam
# do
# 	f2=$(basename $f .hg38_core.merged_filtered.bam)
# 	samtools mpileup --ff UNMAP -B -Q 0 -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/mpileup_out/WGS/${f2}.mpileup_out.txt &
# done
# wait


# for f in $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/*.mpileup_out.txt
# do
# 	gzip -f $f &
# done

# for f in $DATA_PATH/bowtie2_out/mpileup_out/WGS/*.mpileup_out.txt
# do
# 	gzip -f $f &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/*.mpileup_out.txt.gz
# do
# 	f2=$(basename $f .mpileup_out.txt.gz)
# 	zcat $f | cut -f 1-5 | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/${f2}.mpileup_filtered.txt &
# done

# for f in $DATA_PATH/bowtie2_out/mpileup_out/WGS/*.mpileup_out.txt.gz
# do
# 	f2=$(basename $f .mpileup_out.txt.gz)
# 	zcat $f | cut -f 1-5 | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/mpileup_out/WGS/${f2}.mpileup_filtered.txt &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/*.mpileup_filtered.txt
# do
# 	gzip -f $f &
# done

# for f in $DATA_PATH/bowtie2_out/mpileup_out/WGS/*.mpileup_filtered.txt
# do
# 	gzip -f $f &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/baseCt/merged/ATACseq
# for f in $DATA_PATH/bowtie2_out/mpileup_out/ATACseq/*.mpileup_filtered.txt.gz
# do
# 	f2=$(basename $f .mpileup_filtered.txt.gz)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/baseCt/merged/ATACseq/$f2 &
# done

# mkdir -p $DATA_PATH/bowtie2_out/baseCt/merged/WGS
# for f in $DATA_PATH/bowtie2_out/mpileup_out/WGS/*.mpileup_filtered.txt.gz
# do
# 	f2=$(basename $f .mpileup_filtered.txt.gz)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/baseCt/merged/WGS/$f2 &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/baseCt/merged/ATACseq/*.baseCt.txt
# do
# 	gzip -f $f &
# done

# for f in $DATA_PATH/bowtie2_out/baseCt/merged/WGS/*.baseCt.txt
# do
# 	gzip -f $f &
# done
# wait







# ## bwa iterative mapping 
# mkdir -p $DATA_PATH/bwa_out/sam
# mkdir -p $DATA_PATH/bwa_out/bam/mapped_only
# mkdir -p $DATA_PATH/bwa_out/bam/merged
# mkdir -p $DATA_PATH/bwa_out/fastq
# # for f in $DATA_PATH/noAdapter/*_R1.noAdapter.fastq.gz
# # do
# # 	f2=$(basename $f _R1.noAdapter.fastq.gz)
# # 	in_fq1=${f2}_R1.noAdapter.fastq.gz
# # 	in_fq2=${f2}_R2.noAdapter.fastq.gz
# # 	sam_prefix=${f2}.hg38_core
# # 	bam_al=${f2}.hg38_core.al.bam
# # 	out_fq=${f2}.hg38_core.un.fastq
# # 	# bwa mem -t $NUM_THREADS -p $BWA_INDEX_PATH/hg38_core $DATA_PATH/noAdapter/${in_fq1} $DATA_PATH/noAdapter/${in_fq2} > $DATA_PATH/bwa_out/${f2}.hg38_core.sam
# # 	# samtools view -@ $(($NUM_THREADS / 2)) -bS $DATA_PATH/bwa_out/sam/${sam_prefix}.sam | samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bwa_out/bam/${sam_prefix}.bam
# # 	# samtools fastq -@ $NUM_THREADS -f 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/fastq/${out_fq}
# # 	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/bam/mapped_only/${bam_al}

# # 	in_fq=${f2}.hg38_core.un.fastq
# # 	sam_prefix=${f2}.un.B3.hg38_core
# # 	bam_al=${f2}.hg38_core.al2.bam
# # 	out_fq=${f2}.hg38_core.un2.fastq
# # 	bwa mem -t $NUM_THREADS -B 3 -p $DATA_PATH/bwa_index/hg38_core $DATA_PATH/bwa_out/fastq/${in_fq} > $DATA_PATH/bwa_out/sam/${sam_prefix}.sam
# # 	samtools view -@ $(($NUM_THREADS / 2)) -bS $DATA_PATH/bwa_out/sam/${sam_prefix}.sam| samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bwa_out/bam/${sam_prefix}.bam
# # 	samtools fastq -@ $NUM_THREADS -f 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/fastq/${out_fq}
# # 	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/bam/mapped_only/${bam_al}

# # 	in_fq=${f2}.hg38_core.un2.fastq
# # 	sam_prefix=${f2}.un2.B2.hg38_core
# # 	bam_al=${f2}.hg38_core.al3.bam
# # 	out_fq=${f2}.hg38_core.un3.fastq
# # 	bwa mem -t $NUM_THREADS -B 2 -p $DATA_PATH/bwa_index/hg38_core $DATA_PATH/bwa_out/fastq/${in_fq} > $DATA_PATH/bwa_out/sam/${sam_prefix}.sam
# # 	samtools view -@ $(($NUM_THREADS / 2)) -bS $DATA_PATH/bwa_out/sam/${sam_prefix}.sam| samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bwa_out/bam/${sam_prefix}.bam
# # 	samtools fastq -@ $NUM_THREADS -f 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/fastq/${out_fq}
# # 	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/bam/mapped_only/${bam_al}

# # 	in_fq=${f2}.hg38_core.un3.fastq
# # 	sam_prefix=${f2}.un3.B1.hg38_core
# # 	bam_al=${f2}.hg38_core.al4.bam
# # 	out_fq=${f2}.hg38_core.un4.fastq
# # 	bwa mem -t $NUM_THREADS -B 1 -p $DATA_PATH/bwa_index/hg38_core $DATA_PATH/bwa_out/fastq/${in_fq} > $DATA_PATH/bwa_out/sam/${sam_prefix}.sam
# # 	samtools view -@ $(($NUM_THREADS / 2)) -bS $DATA_PATH/bwa_out/sam/${sam_prefix}.sam| samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bwa_out/bam/${sam_prefix}.bam
# # 	samtools fastq -@ $NUM_THREADS -f 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/fastq/${out_fq}
# # 	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/bam/mapped_only/${bam_al}

# # 	in_fq=${f2}.hg38_core.un4.fastq
# # 	sam_prefix=${f2}.un4.B0.hg38_core
# # 	bam_al=${f2}.hg38_core.al5.bam
# # 	out_fq=${f2}.hg38_core.un5.fastq
# # 	bwa mem -t $NUM_THREADS -B 0 -p $DATA_PATH/bwa_index/hg38_core $DATA_PATH/bwa_out/fastq/${in_fq} > $DATA_PATH/bwa_out/sam/${sam_prefix}.sam
# # 	samtools view -@ $(($NUM_THREADS / 2)) -bS $DATA_PATH/bwa_out/sam/${sam_prefix}.sam| samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bwa_out/bam/${sam_prefix}.bam
# # 	samtools fastq -@ $NUM_THREADS -f 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/fastq/${out_fq}
# # 	samtools view -@ $NUM_THREADS -b -F 4 $DATA_PATH/bwa_out/bam/${sam_prefix}.bam > $DATA_PATH/bwa_out/bam/mapped_only/${bam_al}

# # 	samtools merge -@ $NUM_THREADS -f $DATA_PATH/bwa_out/bam/merged/${f2}.hg38_core.merged.bam $DATA_PATH/bwa_out/bam/mapped_only/${f2}.hg38_core.al.bam $DATA_PATH/bwa_out/bam/mapped_only/${f2}.hg38_core.al2.bam $DATA_PATH/bwa_out/bam/mapped_only/${f2}.hg38_core.al3.bam $DATA_PATH/bwa_out/bam/mapped_only/${f2}.hg38_core.al4.bam $DATA_PATH/bwa_out/bam/mapped_only/${f2}.hg38_core.al5.bam
# # 	samtools view -@ $NUM_THREADS -b -f 2 -F 4 $DATA_PATH/bwa_out/bam/merged/${f2}.hg38_core.merged.bam > $DATA_PATH/bwa_out/bam/merged/${f2}.hg38_core.merged_filtered.bam
# # done

# # mkdir -p $DATA_PATH/bwa_out/mpileup_out/
# # for f in $DATA_PATH/bwa_out/bam/merged/*.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bwa_out/mpileup_out/${f2}.iterative.mpileup_out.txt &
# # done
# # wait

# # for f in $DATA_PATH/bwa_out/mpileup_out/*.iterative.mpileup_out.txt.gz
# # do
# # 	f2=$(basename $f .mpileup_out.txt.gz)
# # 	zcat $f | cut -f 1-5 | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bwa_out/mpileup_out/${f2}.mpileup_filtered.txt &
# # done
# # wait

# # mkdir -p $DATA_PATH/bwa_out/baseCt
# # for f in $DATA_PATH/bwa_out/mpileup_out/*.iterative.mpileup_filtered.txt.gz
# # do
# # 	f2=$(basename $f .mpileup_filtered.txt.gz)
# # 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bwa_out/baseCt/$f2 &
# # done
# # wait

# # for f in $DATA_PATH/bwa_out/baseCt/*.baseCt.txt
# # do
# # 	gzip -f $f &
# # done
# # wait


# for f in $DATA_PATH/bwa_out/bam/mapped_only/*.bam
# do
# 	f2=$(basename $f .bam)
# 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bwa_out/mpileup_out/${f2}.iterative.mpileup_out.txt &
# done
# wait

# for f in $DATA_PATH/bwa_out/mpileup_out/*al*.iterative.mpileup_out.txt
# do
# 	f2=$(basename $f .mpileup_out.txt)
# 	cat $f | cut -f 1-5 | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bwa_out/mpileup_out/${f2}.mpileup_filtered.txt &
# done
# wait

# for f in $DATA_PATH/bwa_out/mpileup_out/*.mpileup_filtered.txt
# do
# 	gzip -f $f &
# done
# wait

# mkdir -p $DATA_PATH/bwa_out/baseCt
# for f in $DATA_PATH/bwa_out/mpileup_out/*al*.iterative.mpileup_filtered.txt.gz
# do
# 	f2=$(basename $f .mpileup_filtered.txt.gz)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bwa_out/baseCt/$f2 &
# done
# wait

# for f in $DATA_PATH/bwa_out/baseCt/*.baseCt.txt
# do
# 	gzip -f $f &
# done
# wait













# mkdir -p $DATA_PATH/bowtie2_out/mpileup_out/
# for f in $DATA_PATH/bowtie2_out/bam/mapped_only/*.bam
# do
# 	f2=$(basename $f .bam)
# 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/mpileup_out/${f2}.iterative.mpileup_out.txt &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/mpileup_out/*al*.iterative.mpileup_out.txt
# do
# 	f2=$(basename $f .mpileup_out.txt)
# 	cat $f | cut -f 1-5 | grep -v '+\|-\|*' | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/mpileup_out/${f2}.mpileup_filtered.txt &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/baseCt
# for f in $DATA_PATH/bowtie2_out/mpileup_out/*al*.iterative.mpileup_filtered.txt
# do
# 	f2=$(basename $f .mpileup_filtered.txt)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/baseCt/$f2 &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/baseCt/*al*.baseCt.txt
# do
# 	gzip -f $f &
# done
# wait





# for f in $DATA_PATH/bowtie2_out/070323_0uM_Ddda_ACCESS.hg38_core.sam
# do
# 	f2=$(basename $f .sam)
# 	samtools view -@ $(($NUM_THREADS / 2)) -bS $f | samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bowtie2_out/${f2}.bam
# done

# mkdir -p $DATA_PATH/bowtie2_out/mpileup_out
# for f in $DATA_PATH/bowtie2_out/070323_0uM_Ddda_ACCESS.hg38_core.bam
# do
# 	f2=$(basename $f .hg38_core.bam)
# 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/mpileup_out/${f2}.mpileup_out.txt &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/mpileup_out/070323_0uM_Ddda_ACCESS.mpileup_out.txt
# do
# 	f2=$(basename $f .mpileup_out.txt)
# 	cat $f | cut -f 1-5 | grep -v '+\|-\|*' | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/mpileup_out/${f2}.mpileup_filtered.txt &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/baseCt
# for f in $DATA_PATH/bowtie2_out/mpileup_out/070323_0uM_Ddda_ACCESS.mpileup_filtered.txt
# do
# 	f2=$(basename $f .mpileup_filtered.txt)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/baseCt/$f2 &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/baseCt/*.baseCt.txt
# do
# 	gzip $f &
# done
# wait



# ## remap pipeline for Bowtie2
# mkdir -p $DATA_PATH/bowtie2_out/remap
# for f in $DATA_PATH/bowtie2_out/fastq/*.un.fastq
# do
# 	f2=$(basename $f .fastq)
# 	echo $f2
# 	bowtie2 -p $NUM_THREADS --very-sensitive --score-min L,-0.6,-0.7 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.L07.hg38_core.sam
# 	bowtie2 -p $NUM_THREADS --very-sensitive --score-min L,-0.6,-0.9 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.L09.hg38_core.sam
# 	bowtie2 -p $NUM_THREADS --very-sensitive --score-min L,-0.6,-1.1 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.L11.hg38_core.sam
# 	bowtie2 -p $NUM_THREADS --very-sensitive --mp 5 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.mp5.hg38_core.sam
# 	bowtie2 -p $NUM_THREADS --very-sensitive --mp 3 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.mp3.hg38_core.sam
# 	bowtie2 -p $NUM_THREADS --very-sensitive --mp 1 -x $BOWTIE2_INDEX_PATH/hg38_core -U $f -S $DATA_PATH/bowtie2_out/remap/${f2}.mp1.hg38_core.sam
# done

# for f in $DATA_PATH/bowtie2_out/remap/*.sam
# do
# 	f2=$(basename $f .sam)
# 	samtools view -@ $(($NUM_THREADS / 2)) -bS $f | samtools sort -@ $(($NUM_THREADS / 2)) - -o $DATA_PATH/bowtie2_out/remap/${f2}.bam
# done

# mkdir -p $DATA_PATH/bowtie2_out/remap/mpileup_out
# for f in $DATA_PATH/bowtie2_out/remap/*.bam
# do
# 	f2=$(basename $f .hg38_core.bam)
# 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/remap/mpileup_out/${f2}.mpileup_out.txt &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/remap/mpileup_out/*.mpileup_out.txt
# do
# 	f2=$(basename $f .mpileup_out.txt)
# 	cat $f | cut -f 1-5 | grep -v '+\|-\|*' | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/remap/mpileup_out/${f2}.mpileup_filtered.txt &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/remap/baseCt
# for f in $DATA_PATH/bowtie2_out/remap/mpileup_out/*.mpileup_filtered.txt
# do
# 	f2=$(basename $f .mpileup_filtered.txt)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/remap/baseCt/$f2 &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/remap/baseCt/*.baseCt.txt
# do
# 	gzip $f &
# done
# wait




# ## remap pipeline for BWA
# # mkdir -p $DATA_PATH/bwa_index
# # bwa index -p $DATA_PATH/bwa_index/hg38_core $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa

# mkdir -p $DATA_PATH/bowtie2_out/remap
# for f in $DATA_PATH/bowtie2_out/fastq/*.un.fastq
# do
# 	f2=$(basename $f .fastq)
# 	echo $f2
# 	bwa mem -t 20 -T 20 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.T20.hg38_core.sam
# 	bwa mem -t 20 -T 10 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.T10.hg38_core.sam
# 	bwa mem -t 20 -T 0 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.T0.hg38_core.sam
# 	bwa mem -t 20 -B 3 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.B3.hg38_core.sam
# 	bwa mem -t 20 -B 2 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.B2.hg38_core.sam
# 	bwa mem -t 20 -B 1 -p $DATA_PATH/bwa_index/hg38_core $f > $DATA_PATH/bowtie2_out/remap/${f2}.B1.hg38_core.sam
# done

# for f in $DATA_PATH/bowtie2_out/remap/*.sam
# do
# 	f2=$(basename $f .sam)
# 	samtools view -@ 10 -bS $f | samtools sort -@ 10 - -o $DATA_PATH/bowtie2_out/remap/${f2}.bam
# done

# mkdir -p $DATA_PATH/bowtie2_out/remap/mpileup_out
# for f in $DATA_PATH/bowtie2_out/remap/*.bam
# do
# 	f2=$(basename $f .hg38_core.bam)
# 	samtools mpileup -f $PROJ_PATH/data/ATACseq_041023/genome/hg38_core.fa -d 100000 $f > $DATA_PATH/bowtie2_out/remap/mpileup_out/${f2}.mpileup_out.txt &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/remap/mpileup_out/*.mpileup_out.txt
# do
# 	f2=$(basename $f .mpileup_out.txt)
# 	cat $f | cut -f 1-5 | grep -v '+\|-\|*' | grep -P '\d+\t[cCgG]\t\d+' > $DATA_PATH/bowtie2_out/remap/mpileup_out/${f2}.mpileup_filtered.txt &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/remap/baseCt
# for f in $DATA_PATH/bowtie2_out/remap/mpileup_out/*.mpileup_filtered.txt
# do
# 	f2=$(basename $f .mpileup_filtered.txt)
# 	perl $SCRIPT_PATH/ATACseq_mpileup2baseCt_byChr.pl $f $DATA_PATH/bowtie2_out/remap/baseCt/$f2 &
# done
# wait

# for f in $DATA_PATH/bowtie2_out/remap/baseCt/*.baseCt.txt
# do
# 	gzip $f &
# done
# wait

# mkdir -p $DATA_PATH/bowtie2_out/remap/fastq
# for f in $DATA_PATH/bowtie2_out/remap/*.bam
# do
# 	f2=$(basename $f .bam)
# 	samtools fastq -@ 20 -f 4 $f > $DATA_PATH/bowtie2_out/remap/fastq/${f2}.un.fastq
# 	samtools fastq -@ 20 -F 4 $f > $DATA_PATH/bowtie2_out/remap/fastq/${f2}.al.fastq
# done

# mkdir -p $DATA_PATH/bowtie2_out/remap/fasta
# for f in $DATA_PATH/bowtie2_out/remap/fastq/*.fastq
# do
# 	f2=$(basename $f .fastq)
# 	perl $SCRIPT_PATH_2/fq2fa_collapse_prefix.pl $f $DATA_PATH/bowtie2_out/remap/fasta/${f2}.fa $f2 &
# done
# wait





