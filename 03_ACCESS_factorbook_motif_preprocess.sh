#!/bin/bash -l

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
SCRIPT_PATH_3=/net/bgm/sherwood/UltimaGen_scripts
qsub_options="-N factorbook_preprocess -b y -j y -q all.q@compute-0-0,all.q@compute-0-1 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out"


DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged/all/best_quality
FACTORBOOK_PATH=/net/bgm/sherwood/factorbook_data
ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"

mkdir -p $FACTORBOOK_PATH/UltimaGen/K562
mkdir -p $FACTORBOOK_PATH/UltimaGen/HepG2
mkdir -p $FACTORBOOK_PATH/UltimaGen/stats_part

# step 1
step=20
for i in {0..103}
do
	echo $(($step*$i)) - $(($step*($i+1)-1))
	qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/preprocess_motif_beds_step1.py $FACTORBOOK_PATH/ $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $FACTORBOOK_PATH/UltimaGen/ --start_idx $(($step*$i)) --end_idx $(($step*($i+1)-1)) --stats_out $FACTORBOOK_PATH/UltimaGen/stats_part/preprocess_stats.pt${i}.csv"
done

# step 2
mkdir -p $FACTORBOOK_PATH/UltimaGen/stats
qsub_options="-N factorbook_preprocess_step2 -b y -j y -q all.q@compute-0-0,all.q@compute-0-1 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out"
qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/preprocess_motif_beds_step2.py $FACTORBOOK_PATH/factorbook_MEME_motif_summary_with_chipseq_peaks.csv $FACTORBOOK_PATH/UltimaGen/ $FACTORBOOK_PATH/UltimaGen/stats_part/preprocess_stats $FACTORBOOK_PATH/UltimaGen/stats/TF_summary 0 103 --motif_list_in $FACTORBOOK_PATH/TF_master_list.csv"






