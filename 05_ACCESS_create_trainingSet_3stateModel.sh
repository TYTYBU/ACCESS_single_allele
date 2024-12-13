#!/bin/bash -l

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
SCRIPT_PATH_3=/net/bgm/sherwood/UltimaGen_scripts
qsub_options="-N create_trainingSet_3stateModel -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 16"


DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged/all/best_quality
GENOME_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/genome
FACTORBOOK_PATH=/net/bgm/sherwood/factorbook_data
BED_PATH=$FACTORBOOK_PATH/UltimaGen
H5_PATH=$DATA_PATH/trainingSet_3stateModel_v2
MOTIF_FEATURE_PATH=$FACTORBOOK_PATH/UltimaGen/motif_features_v2
PREDEFINED_LABEL_PATH=$FACTORBOOK_PATH/UltimaGen/predefined_labels_v2
ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"



TF_summary=$FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered.csv
mkdir -p $H5_PATH
mkdir -p $MOTIF_FEATURE_PATH
mkdir -p $PREDEFINED_LABEL_PATH

[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    if [ "$cell_type" == "K562" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/create_trainingSet_3stateModel_v2.py $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURE_PATH $cell_type $TF_name $motif_str $H5_PATH --predefined_labels_outDir $PREDEFINED_LABEL_PATH"
    fi

    if [ "$cell_type" == "HepG2" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
    	qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/create_trainingSet_3stateModel_v2.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $MOTIF_FEATURE_PATH $cell_type $TF_name $motif_str $H5_PATH --predefined_labels_outDir $PREDEFINED_LABEL_PATH"
    fi
done

















