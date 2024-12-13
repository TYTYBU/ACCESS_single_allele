import os
import argparse
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics import precision_recall_curve, average_precision_score, roc_curve, auc, f1_score
import tensorflow as tf
from access_util import *
from plotnine import *
import logging
import pickle

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='This script trains 3-state model to predict unbound/bound/recently bound reads.',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('bamFile', 
                        type=str, 
                        default=None,
                        help=('BAM file containing the aligned reads. \n'
                              'Default: None'))

    parser.add_argument('fastaFile', 
                        type=str, 
                        default=None,
                        help=('Indexed fasta file containing hg38 genome. \n'
                              'Default: None'))

    parser.add_argument('featuresDir', 
                        type=str, 
                        default=None,
                        help=('Motif feature part files directory. \n'
                              'Default: None'))

    parser.add_argument('h5Dir', 
                        type=str, 
                        default=None,
                        help=('Path to h5ad directory containing training and testing data. \n'
                              'Default: None'))

    parser.add_argument('modelDir',
                        type=str,
                        default=None,
                        help=('Output directory for individual model dictionary in pickle format. \n'
                              'Default: None'))

    parser.add_argument('cell_type', 
                        type=str, 
                        default=None,
                        help=('TF cell type. \n'
                              'Default: None'))

    parser.add_argument('TF_name', 
                        type=str, 
                        default=None,
                        help=('TF name. \n'
                              'Default: None'))

    parser.add_argument('motif_str', 
                        type=str, 
                        default=None,
                        help=('TF motif string. \n'
                              'Default: None'))

    parser.add_argument('evaluation_outDir',
                        type=str,
                        default=None,
                        help=('Output directory for evaulation stat files. \n'
                              'Default: None'))

    parser.add_argument('--plot_outDir',
                        type=str,
                        default=None,
                        help=('Output directory for evaluation plots. \n'
                              'Default: None'))

    parser.add_argument('--flank_len',
                        type=int,
                        default=150,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 150'))

    parser.add_argument('--motif_ct_thres_by_bin',
                        type=int,
                        default=1000,
                        help=('Max number of reads to for each chromatin accessibility bin. \n'
                              'Default: 1000'))

    parser.add_argument('--read_ct_thres',
                        type=int,
                        default=50,
                        help=('Minimum number of mapped reads for a motif to be eligible to evaluation. \n'
                              'Default: 50'))
    
    parser.add_argument('--full_model_only',
                        action='store_true',
                        default=False,
                        help='Only train full model. \n'
                             'Default: False')

    return parser.parse_args()

# load arguments
args = parse_args()
bam_path = args.bamFile
fasta_path=args.fastaFile
featuresDir = args.featuresDir
h5Dir = args.h5Dir
modelDir = args.modelDir
cell_type = args.cell_type
TF_name = args.TF_name
motif_str = args.motif_str
evaluation_outDir = args.evaluation_outDir
plot_outDir = args.plot_outDir
flank_len = args.flank_len
motif_ct_thres_by_bin = args.motif_ct_thres_by_bin
read_ct_thres = args.read_ct_thres


# load bam file
logging.info(f'Loading bam file from {bam_path}')
if not os.path.exists(bam_path + '.bai'):
    logging.info(f'Creating index for bam file.')
    pysam.index(bam_path)
bam = pysam.AlignmentFile(bam_path,'rb')

# load genome file
logging.info(f'Loading genome from {fasta_path}')
if not os.path.exists(fasta_path + '.fai'):
    logging.info("fasta file need to indexed first!")
    exit(1)
fasta = pysam.FastaFile(fasta_path)

# load read labeling feature parameters
logging.info(f'Loading motif features from {featuresDir}')
feature_path = featuresDir + f'/{cell_type}_{TF_name}_{motif_str}.motif_features.csv'
df_motif_features = pd.read_csv(feature_path)
df_motif_features['center_pos_list'] = df_motif_features['center_pos_list'].apply(lambda x: ast.literal_eval(x))
df_motif_features['l_peak_pos_list'] = df_motif_features['l_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
df_motif_features['r_peak_pos_list'] = df_motif_features['r_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
df_active_motif_features_expanded = expand_active_motif_features(df_motif_features)

center_trough_pos_list = df_motif_features['center_pos_list'][0]
left_peak_pos_list = df_motif_features['l_peak_pos_list'][0]
right_peak_pos_list = df_motif_features['r_peak_pos_list'][0]

peak_w = df_motif_features['peak_weight'][0]
footprint_w = df_motif_features['footprint_weight'][0]
global_std_w = df_motif_features['global_std_weight'][0]
logging.info(f'Footprint weight: {footprint_w}, Peak weight: {peak_w}, Global std weight: {global_std_w}')

peak_thres = df_motif_features['left_peak'][0]*peak_w*0.5 + df_motif_features['right_peak'][0]*peak_w*0.5 + df_motif_features['footprint'][0]*(1-peak_w)
trough_thres = df_motif_features['left_peak'][0]*footprint_w*0.5 + df_motif_features['right_peak'][0]*footprint_w*0.5 + df_motif_features['footprint'][0]*(1-footprint_w)
global_thres = df_motif_features['inactive_mean'][0] + df_motif_features['inactive_std'][0]*global_std_w


df_prc_all_summary = pd.DataFrame()
df_roc_all_summary = pd.DataFrame()
df_composite_edit_frac = pd.DataFrame()
df_read_type_frac = pd.DataFrame()

model_type_list = ['base_channel', 'edit_channel', 'all_channel']
if args.full_model_only:
    model_type_list = ['all_channel']
read_type_list = ['unbound', 'bound', 'recently bound']

# load from h5py
h5ad_path = h5Dir + f'/{cell_type}_{TF_name}_{motif_str}__3readTypes_oneHot.h5'
logging.info(f'Loading h5ad file from {h5ad_path} ...')
hf = h5py.File(h5ad_path, 'r')
x_test = np.array(hf.get('x_test'))
y_test_class = np.array(hf.get('y_test_class'))
df_clean_motifs = pd.read_hdf(h5ad_path, key='df_clean_motifs')
df_clean_motifs.loc[df_clean_motifs['chipseq_signalVal']<0, 'chipseq_signalVal'] = 0
df_clean_motifs.loc[df_clean_motifs['chipseq_qVal']<0, 'chipseq_qVal'] = 0
df_clean_motifs['chipseq_norm_signalVal'] = df_clean_motifs['chipseq_signalVal'] / df_clean_motifs['chipseq_signalVal'].max()
df_active = df_clean_motifs.loc[(df_clean_motifs['motif_type'] == 'Active') & (df_clean_motifs['use_type'] == 'test') & (df_clean_motifs['motif_read_ct'] >= read_ct_thres)].copy()
df_inactive = df_clean_motifs.loc[(df_clean_motifs['motif_type'] == 'Inactive') & (df_clean_motifs['use_type'] == 'test') & (df_clean_motifs['motif_read_ct'] >= read_ct_thres)].copy()
hf.close()


# add one dimension for conv layer
x_test_conv = np.expand_dims(x_test, axis=3) 

# Categorically encode labels
NUM_CLASSES = np.unique(y_test_class).shape[0]
y_test_oneHot = tf.keras.utils.to_categorical(y_test_class, NUM_CLASSES)

# load models
pkl_path = modelDir + f'/{cell_type}_{TF_name}_{motif_str}__3stateModels.pkl'
logging.info(f'Loading 3-state model from {pkl_path} ...')
with open(pkl_path, 'rb') as file:
    models_dict = pickle.load(file)

logging.info(f'Calculating PRC, ROC and composite edit fractions ...')

# composite edit fractions based on defined labels
obv_classes = np.argmax(y_test_oneHot, axis=1)
df_temp = oneHotMat_2_edit_frac_by_read_type_labels(x_test, obv_classes)
df_temp['label'] = 'pre_defined'
df_temp['cell_type'] = cell_type
df_temp['TF_name'] = TF_name
df_temp['motif_str'] = motif_str
df_composite_edit_frac = pd.concat([df_composite_edit_frac, df_temp], ignore_index=True)

# evaluate model performance on testing data
for model_type in model_type_list:
    if model_type == 'base_channel':
        x_test2 = x_test_conv[:,0:4,:,:]
    elif model_type == 'edit_channel':
        x_test2 = x_test_conv[:,4:7,:,:]
    elif model_type == 'all_channel':
        x_test2 = x_test_conv

    y_test_pred = models_dict[model_type].predict(x_test2, verbose=0)
    y_pred_class = np.argmax(y_test_pred, axis=1)
    f1 = f1_score(y_test_class, y_pred_class, average='weighted')
    logging.info(f'{model_type} F1 score (weighted): {f1:.3f}')

    # read composite edit fractions based on predicted labels
    df_temp = oneHotMat_2_edit_frac_by_read_type_labels(x_test, y_pred_class)
    df_temp['label'] = model_type
    df_temp['cell_type'] = cell_type
    df_temp['TF_name'] = TF_name
    df_temp['motif_str'] = motif_str
    df_composite_edit_frac = pd.concat([df_composite_edit_frac, df_temp], ignore_index=True)

    for i in range(len(read_type_list)):
        precision, recall, thres = precision_recall_curve(y_test_oneHot[:,i], y_test_pred[:,i])
        # average_precision = average_precision_score(y_test_oneHot[:,i], y_test_pred[:,i])
        df_prc = pd.DataFrame({'precision': precision, 'recall': recall})
        auprc = auc(recall, precision)
        df_prc['read_type'] = read_type_list[i]
        df_prc['AUPRC'] = auprc
        df_prc['f1_weighted'] = f1
        df_prc['model_type'] = model_type
        df_prc['cell_type'] = cell_type
        df_prc['TF_name'] = TF_name
        df_prc['motif_str'] = motif_str
        df_prc_summary = df_prc[['cell_type', 'TF_name', 'motif_str', 'model_type', 'f1_weighted', 'read_type', 'AUPRC']].drop_duplicates(ignore_index=True)
        df_prc_all_summary = pd.concat([df_prc_all_summary, df_prc_summary], ignore_index=True)

        roc_fpr, roc_tpr, roc_thres = roc_curve(y_test_oneHot[:,i], y_test_pred[:,i])
        df_roc = pd.DataFrame({'fpr': roc_fpr, 'tpr': roc_tpr, 'threshold': roc_thres})
        roc_auc = auc(roc_fpr, roc_tpr)
        df_roc['read_type'] = read_type_list[i]
        df_roc['AUROC'] = roc_auc
        df_roc['model_type'] = model_type
        df_roc['cell_type'] = cell_type
        df_roc['TF_name'] = TF_name
        df_roc['motif_str'] = motif_str
        df_roc_summary = df_roc[['cell_type', 'TF_name', 'motif_str', 'model_type', 'read_type', 'AUROC']].drop_duplicates(ignore_index=True)
        df_roc_all_summary = pd.concat([df_roc_all_summary, df_roc_summary], ignore_index=True)


logging.info(f'Calculating read type fractions ...')

# bin test data motifs by chipseq scores
score_lower_bound_list = [x for x in range(0,401,50)]
score_upper_bound_list = score_lower_bound_list[1:] + [float('inf')]
score_bins = ['Inactive']
df_motifs_bin = df_inactive.sample(motif_ct_thres_by_bin) if len(df_inactive)>motif_ct_thres_by_bin else df_inactive
df_motifs_bin['score_bin'] = 'Inactive'
for score_lower_bound, score_upper_bound in zip(score_lower_bound_list, score_upper_bound_list):
    score_bin = 'Active: ' + str(score_lower_bound) + ('+' if score_upper_bound == float('inf') else ('-' + str(score_upper_bound)))
    score_bins.append(score_bin)
    df_temp = df_active.loc[df_active['chipseq_signalVal'].between(score_lower_bound, score_upper_bound, inclusive='left')].copy()
    df_temp['score_bin'] = score_bin
    df_temp = df_temp.sample(motif_ct_thres_by_bin) if len(df_temp)>motif_ct_thres_by_bin else df_temp
    df_motifs_bin = pd.concat([df_motifs_bin, df_temp], ignore_index=True)

df_motifs_ct = df_motifs_bin.groupby('score_bin').size().reset_index()
df_motifs_ct['score_bin'] = pd.Categorical(df_motifs_ct['score_bin'], categories=score_bins)
df_motifs_ct.sort_values('score_bin', inplace=True)
df_motifs_ct.rename(columns={0: 'motif_ct'}, inplace=True)

# calculate read fractions of 3 read types for each motif
for score_bin in score_bins:
    df_motifs_bin_temp = df_motifs_bin.loc[df_motifs_bin['score_bin'] == score_bin].copy()
    logging.info(f'{score_bin} motif count: {len(df_motifs_bin_temp)}')

    for index, row in df_motifs_bin_temp.iterrows():
        oneHotMat_dict = motif_row_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, row, 
                                                                  center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                                                  trough_thres, peak_thres, global_thres,
                                                                  flank_len = flank_len)
        oneHotMat_row = np.array([])
        obv_classes_row = []
        for i in range(len(read_type_list)):
            oneHotMat = oneHotMat_dict[read_type_list[i]]
            obv_classes = [i] * oneHotMat_dict[read_type_list[i]].shape[0]
            if oneHotMat.size > 0:
                if oneHotMat_row.size == 0:
                    oneHotMat_row = oneHotMat
                    obv_classes_row = obv_classes
                else:
                    oneHotMat_row = np.vstack((oneHotMat_row, oneHotMat))
                    obv_classes_row = obv_classes_row + obv_classes
        
        if oneHotMat_row.shape[0] > 0:
            # use pre-defined labels as read type
            uniq_value, value_ct = np.unique(obv_classes_row, return_counts=True)
            df_temp2 = pd.DataFrame({'read_type': [read_type_list[x] for x in uniq_value], 
                                    'read_ct_defined': value_ct, 'read_frac_defined': value_ct/np.sum(value_ct)})

            # use model to predict read type
            oneHotMat_conv = np.expand_dims(oneHotMat_row, axis=3)
            pred_probs = models_dict['all_channel'].predict(oneHotMat_conv, verbose=0)
            pred_classes = np.argmax(pred_probs, axis=1)
            mean_probs = np.mean(pred_probs, axis=0)

            uniq_value, value_ct = np.unique(pred_classes, return_counts=True)
            df_temp3 = pd.DataFrame({'read_type': [read_type_list[x] for x in uniq_value], 
                                    'read_ct_predicted': value_ct, 'read_frac_predicted': value_ct/np.sum(value_ct)})

            df_temp5 = pd.DataFrame({'read_type': read_type_list, 'read_frac_predicted_prob': mean_probs})
            df_temp4 = pd.merge(df_temp2, df_temp3, how='outer', on=['read_type'])
            df_temp4 = pd.merge(df_temp4, df_temp5, how='outer', on='read_type')
            df_temp4['score_bin'] = score_bin
            df_temp4['chipseq_signalVal'] = row['chipseq_signalVal']
            df_temp4['chipseq_norm_signalVal'] = row['chipseq_norm_signalVal']
            df_temp4['cell_type'] = cell_type
            df_temp4['TF_name'] = TF_name
            df_temp4['motif_str'] = motif_str

            df_read_type_frac = pd.concat([df_read_type_frac, df_temp4], ignore_index=True)


if evaluation_outDir is not None:
    logging.info(f'Output evaluation stats to {evaluation_outDir}')

    out_path = evaluation_outDir + f'/{cell_type}_{TF_name}_{motif_str}.prc.csv'
    df_prc_all_summary.to_csv(out_path, index=False)
    out_path = evaluation_outDir + f'/{cell_type}_{TF_name}_{motif_str}.roc.csv'
    df_roc_all_summary.to_csv(out_path, index=False)
    out_path = evaluation_outDir + f'/{cell_type}_{TF_name}_{motif_str}.edit_frac.csv'
    df_composite_edit_frac.to_csv(out_path, index=False)
    out_path = evaluation_outDir + f'/{cell_type}_{TF_name}_{motif_str}.read_type_frac.csv'
    df_read_type_frac.to_csv(out_path, index=False)
    out_path = evaluation_outDir + f'/{cell_type}_{TF_name}_{motif_str}.chromAcc_bin_ct.csv'
    df_motifs_ct.to_csv(out_path, index=False)


if plot_outDir is not None:
    logging.info(f'Output evaluation plots to {plot_outDir}')

    # edit fraction plot
    df_composite_edit_frac['label'] = pd.Categorical(df_composite_edit_frac['label'], categories=['pre_defined', 'base_channel', 'edit_channel', 'all_channel'])
    plt = (
        ggplot(df_composite_edit_frac, aes(x='relative_pos', y='edit_frac', color='read_type'))
        + geom_point(size=0.3, alpha=0.2)
        + stat_smooth(method='mavg', method_args={'window': 5}, size=0.25, se=False)
        + labs(x='Relative position', y='Edit fraction')
        + facet_wrap('label', nrow=1) + ggtitle(f'{cell_type} {TF_name} {motif_str}')
        + theme_light()
    )
    plt_path = plot_outDir + f'/{cell_type}_{TF_name}_{motif_str}.edit_frac.pdf'
    if args.full_model_only:
        plt.save(filename=plt_path, width=6, height=3)
    else:
        plt.save(filename=plt_path, width=12, height=3)

    # violin plot
    df_read_type_frac['read_type'] = pd.Categorical(df_read_type_frac['read_type'], categories=read_type_list)
    df_read_type_frac['score_bin'] = pd.Categorical(df_read_type_frac['score_bin'], categories=score_bins)
    df_read_type_frac_long = df_read_type_frac.melt(id_vars=['cell_type', 'TF_name', 'motif_str', 'read_type', 'score_bin'], 
                                                    value_vars=['read_frac_defined', 'read_frac_predicted', 'read_frac_predicted_prob'], 
                                                    var_name='read_frac_type', value_name='read_frac_value')
    df_read_type_frac_long['read_frac_type'] = df_read_type_frac_long['read_frac_type'].str.replace('read_frac_', '')
    plt = (
        ggplot(df_read_type_frac_long, aes(x='score_bin', y='read_frac_value', fill='score_bin')) +
        geom_violin(size=0.25) + geom_boxplot(fill=None, width=0.2, outlier_alpha=0.25, size=0.25, outlier_size=0.25) + 
        facet_grid('read_frac_type ~ read_type') + theme_light() + ggtitle(f'{cell_type} {TF_name} {motif_str}') +
        theme(axis_text_x=element_text(rotation=45, hjust=1))
    )
    plt_path = plot_outDir + f'/{cell_type}_{TF_name}_{motif_str}.read_type_frac.violinPlt.pdf'
    plt.save(filename=plt_path, width=12, height=8)

    # correlation plot 
    df_read_type_frac['read_type'] = pd.Categorical(df_read_type_frac['read_type'], categories=read_type_list)
    df_read_type_frac2 = df_read_type_frac.loc[df_read_type_frac['chipseq_norm_signalVal'] > 0].copy()
    df_read_type_frac2['chipseq_norm_signalVal_log10'] = np.log10(df_read_type_frac2['chipseq_norm_signalVal'])
    df_read_type_frac2_long = df_read_type_frac2.melt(id_vars=['cell_type', 'TF_name', 'motif_str', 'read_type', 'chipseq_norm_signalVal_log10'], value_vars=['read_frac_defined', 'read_frac_predicted', 'read_frac_predicted_prob'], var_name='read_frac_type', value_name='read_frac_value')
    df_read_type_frac2_long['read_frac_type'] = df_read_type_frac2_long['read_frac_type'].str.replace('read_frac_', '')

    df_spearman_corr = df_read_type_frac2_long.groupby(['TF_name', 'read_frac_type', 'read_type']).apply(lambda grp: pd.Series({**spearmanr_wrapper(grp['read_frac_value'], grp['chipseq_norm_signalVal_log10'], prefix='spearman_'), **pearsonr_wrapper(grp['read_frac_value'], grp['chipseq_norm_signalVal_log10'], prefix='pearson_')}), include_groups=False).reset_index()
    df_spearman_corr['corr_str'] = df_spearman_corr.apply(lambda row: f'rho={row['spearman_rho']:.2f}, p={row['spearman_pval']:.2e}\nr={row['pearson_r']:.2f}, p={row['pearson_pval']:.2e}', axis=1)
    df_temp = df_read_type_frac2_long.groupby(['read_frac_type', 'read_type']).apply(lambda grp: pd.Series({'x': grp['read_frac_value'].max(), 'y': grp['chipseq_norm_signalVal_log10'].max()}), include_groups=False).reset_index()
    df_spearman_corr = pd.merge(df_spearman_corr, df_temp, on=['read_frac_type', 'read_type'])

    plt = (
        ggplot() + 
        geom_point(df_read_type_frac2_long, aes(x='read_frac_value', y='chipseq_norm_signalVal_log10', color='chipseq_norm_signalVal_log10'), alpha=0.25, size=0.25) + 
        scale_color_gradient(low="blue", high="red") + 
        geom_text(df_spearman_corr, aes(label='corr_str', x='x', y='y'), ha='right', va='top') + 
        theme_light() + facet_grid('read_frac_type ~ read_type', scales='free') + labs(color='chipseq_norm_signalVal_log10') + 
        ggtitle(f'{cell_type} {TF_name} {motif_str}') +
        theme(legend_position='bottom')
    )
    plt_path = plot_outDir + f'/{cell_type}_{TF_name}_{motif_str}.read_type_frac.corPlt.pdf'
    plt.save(filename=plt_path, width=9, height=9.5)


    


