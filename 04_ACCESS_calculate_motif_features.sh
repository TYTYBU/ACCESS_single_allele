#!/bin/bash -l

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
SCRIPT_PATH_3=/net/bgm/sherwood/UltimaGen_scripts


DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged/all/best_quality
GENOME_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/genome
FACTORBOOK_PATH=/net/bgm/sherwood/factorbook_data
ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"


## step 1: calculate motif features
qsub_options="-N calc_motif_features -b y -j y -q all.q@compute-0-0,all.q@compute-0-1 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 1"

TF_summary=$FACTORBOOK_PATH/UltimaGen/stats/TF_summary.selected.csv
MOTIF_FEATURES_PATH=$FACTORBOOK_PATH/UltimaGen/motif_features_v2
mkdir -p $MOTIF_FEATURES_PATH

[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    if [ "$cell_type" == "K562" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/calc_motif_features_step1_v2.py $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $FACTORBOOK_PATH/UltimaGen $cell_type $TF_name $motif_str $MOTIF_FEATURES_PATH --sampled_motif_ct 2000"
    fi

    if [ "$cell_type" == "HepG2" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/calc_motif_features_step1_v2.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $FACTORBOOK_PATH/UltimaGen $cell_type $TF_name $motif_str $MOTIF_FEATURES_PATH --sampled_motif_ct 2000"
    fi
done

## step 2: combined all motif features
qsub_options="-N calc_motif_features_step2 -b y -j y -q all.q@compute-0-0,all.q@compute-0-1 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out"
qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/calc_motif_features_step2.py $FACTORBOOK_PATH/UltimaGen/stats/TF_summary.selected.csv $MOTIF_FEATURES_PATH $FACTORBOOK_PATH/UltimaGen/stats"

