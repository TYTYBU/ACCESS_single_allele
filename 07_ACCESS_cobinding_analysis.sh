#!/bin/bash -l

PROJ_PATH=/net/bgm/sherwood/NGS_analysis_proj
SCRIPT_PATH=$PROJ_PATH/script
SCRIPT_PATH_2=/net/home/tianyu/Tian_script
SCRIPT_PATH_3=/net/bgm/sherwood/UltimaGen_scripts


qsub_options="-N cobinding_definedlabel -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 1"
cell_type='HepG2'
TF_name='HLF'
motif_str='ENCSR528PSI_GTTATRCAACH'
qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/cobinding_definedLabel.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $TF_summary $FACTORBOOK_PATH/UltimaGen/ $cell_type $TF_name $motif_str $CODBIND_PATH --read_ct_thres 50 --correlation_thres 0 --max_dist 100"  


#### cobinding analysis with modeled data
DATA_PATH=/net/bgm/sherwood/UltimaGen_datasets/bulkATACseq
BAM_PATH=$DATA_PATH/bowtie2_out/bam/merged/all/best_quality
GENOME_PATH=/net/bgm/sherwood/NGS_analysis_proj/data/ATACseq_041023/genome
FACTORBOOK_PATH=/net/bgm/sherwood/factorbook_data
H5_PATH=$DATA_PATH/trainingSet_3stateModel_v2
MODEL_PATH=$FACTORBOOK_PATH/UltimaGen/model_3stateModel_v2
ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"
qsub_options="-N cobinding_3stateModel_v2 -b y -j y -q all.q@compute-0-0,all.q@compute-0-1,all.q@compute-0-2 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe serial 1"

TF_summary=$FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered_HepG2_modeled.csv
read_ct_thres=50
CODBIND_PATH=$FACTORBOOK_PATH/UltimaGen/cobinding_v2_read_thres50_rho0/HepG2_footprint_mean
mkdir -p $CODBIND_PATH

[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    if [ "$cell_type" == "HepG2" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/cobinding_3stateModel_v2.py $BAM_PATH/062424_HepG2_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $TF_summary $FACTORBOOK_PATH/UltimaGen/ $MODEL_PATH $cell_type $TF_name $motif_str $CODBIND_PATH --read_ct_thres ${read_ct_thres} --correlation_thres 0 --min_dist 10 --max_dist 100 --adj_method_num 3"  
    fi
done

TF_summary=$FACTORBOOK_PATH/UltimaGen/stats/TF_summary.features_filtered_K562_modeled.csv
read_ct_thres=50
CODBIND_PATH=$FACTORBOOK_PATH/UltimaGen/cobinding_v2_read_thres50_rho0/K562_footprint_mean
mkdir -p $CODBIND_PATH

[ ! -f $TF_summary ] && { echo "$TF_summary file not found"; exit 99; }
tail -n +2 "$TF_summary" | while IFS=',' read -r cell_type TF_name motif_str rest
do
    if [ "$cell_type" == "K562" ]; then
        echo "cell type: $cell_type, TF name: $TF_name, motif string: $motif_str"
        qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/cobinding_3stateModel_v2.py $BAM_PATH/062424_K562_concurrent_ACCESS-ATAC.dedup.hg38_core.bam $GENOME_PATH/hg38_core.fa $TF_summary $FACTORBOOK_PATH/UltimaGen/ $MODEL_PATH $cell_type $TF_name $motif_str $CODBIND_PATH --read_ct_thres ${read_ct_thres} --correlation_thres 0 --min_dist 10 --max_dist 100 --adj_method_num 3"  
    fi
done


# ### final co-bidning coditions test
# ENV_SETTING="eval \"\$(/net/apps/conda/miniconda3/bin/conda shell.bash hook)\" && conda activate access"
# qsub_options="-N cobinding_3stateModel_v2 -b y -j y -q all.q@compute-0-1 -wd /net/bgm/sherwood -o /net/bgm/sherwood/job_out -pe smp 16"
# qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/test.py"
# qsub $qsub_options "$ENV_SETTING && python $SCRIPT_PATH/test_2.py"








