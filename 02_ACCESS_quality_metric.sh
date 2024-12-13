#!/bin/bash -l

#$ -N quality_metric
#$ -j y
#$ -V
#$ -q all.q@compute-0-0,all.q@compute-0-1
#$ -pe serial 32
#$ -m ae
#$ -M tyu7@bwh.harvard.edu
#$ -wd /net/bgm/sherwood/NGS_analysis_proj
#$ -o /net/bgm/sherwood/NGS_analysis_proj/job_out

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
GENOME_PATH=$PROJ_PATH/data/ATACseq_041023/genome


DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/hisat3n_out/bam/all_dedup
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.hg38_core.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 180 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/hisat3n_out/bam/all_dedup
# PEAK_PATH=$PROJ_PATH/data/ACCESS_central/factorbook_peaks_default
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.bam
# do
# 	f2=$(basename $f .hg38_core.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate_v2.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_NRF_baseCt.csv --no_ref_n --use_strand &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/hisat3n_out/bam/all_dedup
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.bam
# do
# 	f2=$(basename $f .hg38_core.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# 	wait
# done
# wait


### twin reads only!
# quality metric - TSS enrichment
conda activate tss_enrich
BAM_PATH=$DATA_PATH/hisat3n_out/bam/full_len_dedup/twin_merged
mkdir -p $DATA_PATH/quality_metric/twin_merged/tss_enrich
for f in $BAM_PATH/*.hg38_core.bam
do
	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 180 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/twin_merged/tss_enrich &
done
wait

# quanlity metric - CTCF peak/trough peak/background
conda activate bio
BAM_PATH=$DATA_PATH/hisat3n_out/bam/full_len_dedup/twin_merged
PEAK_PATH=$PROJ_PATH/data/ACCESS_central/factorbook_peaks_default
mkdir -p $DATA_PATH/quality_metric/twin_merged/baseCt
for f in $BAM_PATH/*.hg38_core.bam
do
	f2=$(basename $f .hg38_core.bam)
	python $SCRIPT_PATH/get_motif_edit_rate_v2.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/twin_merged/${f2}.CTCF_NRF_baseCt.csv --no_ref_n --use_strand &
done
wait

# quanlity metric - accessible chromatin read enrichment
conda activate bio
BAM_PATH=$DATA_PATH/hisat3n_out/bam/full_len_dedup/twin_merged
BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

mkdir -p $DATA_PATH/quality_metric/twin_merged/chromeAcc_bin
for f in $BAM_PATH/*.hg38_core.bam
do
	f2=$(basename $f .hg38_core.bam)
	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/twin_merged/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/twin_merged/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/twin_merged/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
	wait
done
wait




DATA_PATH=$PROJ_PATH/data/Rich_ACCESS_062124

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/ACCESS_central/factorbook_peaks_default
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate_v2.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_NRF_baseCt.csv --no_ref_n --use_strand &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# 	wait
# done
# wait



DATA_PATH=$PROJ_PATH/data/Rich_ACCESS_060724

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/ACCESS_central/factorbook_peaks_default
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_NRF_baseCt.csv --no_ref_n &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# 	wait
# done
# wait


DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_051324

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/ACCESS_central/factorbook_peaks_default
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_NRF_baseCt.csv --no_ref_n &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# 	wait
# done
# wait


DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_050924

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# 	wait
# done
# wait



DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_042724

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done
# wait

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done
# wait

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done
# wait



DATA_PATH=$PROJ_PATH/data/Ellie_scACCESS_121823

# # # quality metric - TSS enrichment
# # conda activate tss_enrich
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # mkdir -p $DATA_PATH/quality_metric/tss_enrich
# # for f in $BAM_PATH/*.merged.bam
# # do
# # 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# # done

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/sc_valid/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_030424

# # # quality metric - TSS enrichment
# # conda activate tss_enrich
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # mkdir -p $DATA_PATH/quality_metric/tss_enrich
# # for f in $BAM_PATH/*.merged.bam
# # do
# # 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# # done

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Jake_ACCESS_030424

# # # quality metric - TSS enrichment
# # conda activate tss_enrich
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # mkdir -p $DATA_PATH/quality_metric/tss_enrich
# # for f in $BAM_PATH/*.merged.bam
# # do
# # 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# # done

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_020824

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # quanlity metric - CTCF peak/trough peak/background
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# mkdir -p $DATA_PATH/quality_metric/baseCt
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done





DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_customTn5_020824

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done





DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_012524

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_bulk_121823

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_112223

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged_filtered.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 150 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged_filtered.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged_filtered.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged_filtered.bam
# do
# 	f2=$(basename $f .hg38_core.merged_filtered.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done



DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_111523

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 75 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 75 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*HCT116*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Jake_ACCESS_111523

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 75 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done




DATA_PATH=$PROJ_PATH/data/Ellie_ACCESS_092523

# # quality metric - TSS enrichment
# conda activate tss_enrich
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# mkdir -p $DATA_PATH/quality_metric/tss_enrich
# for f in $BAM_PATH/*.merged.bam
# do
# 	python $SCRIPT_PATH/encode_task_tss_enrich.py --read-len 75 --nodup-bam $f --chrsz $GENOME_PATH/auxiliary_data/hg38.chrom.sizes --tss $GENOME_PATH/auxiliary_data/hg38.refGene.TSS.bed --out-dir $DATA_PATH/quality_metric/tss_enrich &
# done

# # # quanlity metric - CTCF peak/trough peak/background
# # conda activate bio
# # BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# # PEAK_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/CTCF_peaks
# # mkdir -p $DATA_PATH/quality_metric/baseCt
# # for f in $BAM_PATH/*.hg38_core.merged.bam
# # do
# # 	f2=$(basename $f .hg38_core.merged.bam)
# # 	python $SCRIPT_PATH/get_motif_edit_rate.py $f $GENOME_PATH/hg38_core.fa $PEAK_PATH --baseCt_grouped $DATA_PATH/quality_metric/baseCt/${f2}.CTCF_baseCt.csv --no_ref_n &
# # done

# # quanlity metric - accessible chromatin read enrichment
# conda activate bio
# BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged
# BW_PATH=$PROJ_PATH/data/Ellie_ACCESS_102323/Encode_accessibility_scores

# mkdir -p $DATA_PATH/quality_metric/chromeAcc_bin
# for f in $BAM_PATH/*.hg38_core.merged.bam
# do
# 	f2=$(basename $f .hg38_core.merged.bam)
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/K562_ATACseq_ENCFF600FDO.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.K562.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HepG2_ATACseq_ENCFF262URW.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HepG2.csv --segment_len 100000000 --no_ref_n &
# 	python $SCRIPT_PATH/chromAcc_bin_analysis.py $f $GENOME_PATH/hg38_core.fa $BW_PATH/HCT116_ATACseq_ENCFF962GFP.bigWig --acc_bins 0.1,1 --baseCt_grouped $DATA_PATH/quality_metric/chromeAcc_bin/${f2}.chromeAcc_bin_lite.HCT116.csv --segment_len 100000000 --no_ref_n &
# done



wait






