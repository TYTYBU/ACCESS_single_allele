#!/bin/bash -l

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
SCRIPT_PATH_3=/net/bgm/sherwood/UltimaGen_scripts



DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged/all/best_quality
GENOME_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/genome
FACTORBOOK_PATH=/net/bgm/sherwood/factorbook_data
H5_PATH=$DATA_PATH/trainingSet_3stateModel_v2
MOTIF_FEATURES_PATH=$FACTORBOOK_PATH/UltimaGen/motif_features_v2
MODEL_PATH=$FACTORBOOK_PATH/UltimaGen/model_3stateModel_v2
EVAL_PATH=$FACTORBOOK_PATH/UltimaGen/evaluation_3stateModel_v2
ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"



TF_summary=$FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered.csv

## step 1: train models
qsub_options="-N train_3StateModel -b y -j y -q all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 30"
mkdir -p $MODEL_PATH

[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
    qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/train_3StateModel_single_h5.py $H5_PATH $cell_type $TF_name $motif_str $MODEL_PATH --early_termination --num_epochs 500"
    
    # model_fname="$MODEL_PATH/${cell_type}_${TF_name}_${motif_str}__3stateModels.pkl"
    # if [ "$cell_type" == "HepG2" ]; then
    #     if [ ! -f "$model_fname" ]; then
    #         echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
    #         qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/train_3StateModel_single_h5.py $H5_PATH $cell_type $TF_name $motif_str $MODEL_PATH --early_termination --num_epochs 500 --full_model_only"
    #     fi
    # fi

    # if [ "$cell_type" == "K562" ]; then
    #     if [ ! -f "$model_fname" ]; then
    #         echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
    #         qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/train_3StateModel_single_h5.py $H5_PATH $cell_type $TF_name $motif_str $MODEL_PATH --early_termination --num_epochs 500 --full_model_only"
    #     fi
    # fi
done



## step 2: evaluate models 
read_ct_thres=100

qsub_options="-N evaluate_3StateModel -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 8"
mkdir -p $EVAL_PATH/K562_${read_ct_thres}/plots
mkdir -p $EVAL_PATH/HepG2_${read_ct_thres}/plots
[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    if [ "$cell_type" = "K562" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_single_h5_v2.py $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURES_PATH $H5_PATH $MODEL_PATH $cell_type $TF_name $motif_str $EVAL_PATH/K562_${read_ct_thres} --plot_outDir $EVAL_PATH/K562_${read_ct_thres}/plots --read_ct_thres ${read_ct_thres} --full_model_only"

        # out_fname="$FACTORBOOK_PATH/UltimaGen/evaluation_3stateModel/K562_${read_ct_thres}/${cell_type}_${TF_name}_${motif_str}.roc.csv"
        # model_fname="$FACTORBOOK_PATH/UltimaGen/model_3stateModel/${cell_type}_${TF_name}_${motif_str}__3stateModels.pkl"
        # if [ ! -f "$out_fname" ] && [ -f "$model_fname" ]; then
        #     echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        #     qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_single_h5_v2.py $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURES_PATH $H5_PATH $MODEL_PATH $cell_type $TF_name $motif_str $EVAL_PATH/K562_${read_ct_thres} --plot_outDir $EVAL_PATH/K562_${read_ct_thres}/plots --read_ct_thres ${read_ct_thres} --full_model_only"
        # fi
    fi

    if [ "$cell_type" == "HepG2" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_single_h5_v2.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURES_PATH $H5_PATH $MODEL_PATH $cell_type $TF_name $motif_str $EVAL_PATH/HepG2_${read_ct_thres} --plot_outDir $EVAL_PATH/HepG2_${read_ct_thres}/plots --read_ct_thres ${read_ct_thres} --full_model_only"

        # out_fname="$FACTORBOOK_PATH/UltimaGen/evaluation_3stateModel/HepG2_${read_ct_thres}/${cell_type}_${TF_name}_${motif_str}.roc.csv"
        # model_fname="$FACTORBOOK_PATH/UltimaGen/model_3stateModel/${cell_type}_${TF_name}_${motif_str}__3stateModels.pkl"
        # if [ ! -f "$out_fname" ] && [ -f "$model_fname" ]; then
        #     echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        #     qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_single_h5_v2.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURES_PATH $H5_PATH $MODEL_PATH $cell_type $TF_name $motif_str $EVAL_PATH/HepG2_${read_ct_thres} --plot_outDir $EVAL_PATH/HepG2_${read_ct_thres}/plots --read_ct_thres ${read_ct_thres} --full_model_only"
        # fi
    fi
done


## setp3: summarize evaluation results and filter motif set
read_ct_thres=100

# HepG2
mkdir -p $EVAL_PATH/HepG2_${read_ct_thres}/summary
qsub_options="-N evaluate_3StateModel_summary -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out"
qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_summarize.py $EVAL_PATH/HepG2_${read_ct_thres} --summary_eval_statsDir $EVAL_PATH/HepG2_${read_ct_thres}/summary --plot_outDir $EVAL_PATH/HepG2_${read_ct_thres}/summary --motif_stats_path $FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered.csv --cell_type HepG2 --filtered_motif_stats_path $FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered_HepG2_modeled.csv"

# K562
mkdir -p $EVAL_PATH/K562_${read_ct_thres}/summary
qsub_options="-N evaluate_3StateModel_summary -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out"
qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/evaluate_3StateModel_summarize.py $EVAL_PATH/K562_${read_ct_thres} --summary_eval_statsDir $EVAL_PATH/K562_${read_ct_thres}/summary --plot_outDir $EVAL_PATH/K562_${read_ct_thres}/summary --motif_stats_path $FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered.csv --cell_type K562 --filtered_motif_stats_path $FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered_K562_modeled.csv"











