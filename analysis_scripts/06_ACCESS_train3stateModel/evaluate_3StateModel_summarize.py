import os
import argparse
import numpy as np
import pandas as pd
import h5py
from plotnine import *
from access_util import *
import logging
from statsmodels.stats.multitest import multipletests

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
    parser.add_argument('eval_statsDir',
                        type=str,
                        default=None,
                        help=('Directory with model evaluation output files. \n'
                              'Default: None'))
    
    parser.add_argument('--summary_eval_statsDir',
                        type=str,
                        default=None,
                        help=('Directory to output summary evaluation files. \n'
                              'Default: None'))

    parser.add_argument('--plot_outDir',
                        type=str,
                        default=None,
                        help=('Output directory for evaluation plots. \n'
                              'Default: None'))
    
    parser.add_argument('--motif_stats_path',
                        type=str,
                        default=None,
                        help=('input motif stats file path. \n'
                              'Default: None'))
    
    parser.add_argument('--cell_type', 
                        type=str, 
                        default=None,
                        help=('TF cell type. \n'
                              'Default: None'))

    parser.add_argument('--filtered_motif_stats_path',
                        type=str,
                        default=None,
                        help=('Output filtered motif stats file path. \n'
                              'Default: None'))
    
    parser.add_argument('--full_model_only',
                        action='store_true',
                        default=False,
                        help='Only train full model. \n'
                             'Default: False')

    return parser.parse_args()

# load arguments
args = parse_args()
eval_statsDir = args.eval_statsDir
summary_eval_statsDir=args.summary_eval_statsDir
plot_outDir = args.plot_outDir
motif_stats_path = args.motif_stats_path
cell_type = args.cell_type
filtered_motif_stats_path = args.filtered_motif_stats_path

read_type_list = ['unbound', 'bound', 'recently bound']

# plot model evaluation metrics
stats_file_suffix = '.prc.csv'
df_motif_list = pd.DataFrame([x.replace(stats_file_suffix, '').split('_', maxsplit=2) + [x] for x in os.listdir(eval_statsDir) if stats_file_suffix in x], columns=['cell_type', 'TF_name', 'motif_str', 'file_name'])
df_auprc_all = pd.DataFrame()
for index, row in df_motif_list.iterrows():
    df_auprc = pd.read_csv(eval_statsDir + '/' + row['file_name'])
    df_auprc_all = pd.concat([df_auprc_all, df_auprc])
df_f1_all = df_auprc_all[['cell_type', 'TF_name','motif_str', 'model_type', 'f1_weighted']].drop_duplicates(ignore_index=True)

stats_file_suffix = '.roc.csv'
df_motif_list = pd.DataFrame([x.replace(stats_file_suffix, '').split('_', maxsplit=2) + [x] for x in os.listdir(eval_statsDir) if stats_file_suffix in x], columns=['cell_type', 'TF_name', 'motif_str', 'file_name'])
df_auroc_all = pd.DataFrame()
for index, row in df_motif_list.iterrows():
    df_auroc = pd.read_csv(eval_statsDir + '/' + row['file_name'])
    df_auroc_all = pd.concat([df_auroc_all, df_auroc])

# boxplot for the correlation between read type fractions and normalized chipseq sigVal
stats_file_suffix = '.read_type_frac.csv'
df_motif_list = pd.DataFrame([x.replace(stats_file_suffix, '').split('_', maxsplit=2) + [x] for x in os.listdir(eval_statsDir) if stats_file_suffix in x], columns=['cell_type', 'TF_name', 'motif_str', 'file_name'])

df_read_type_frac_all = pd.DataFrame()
for index, row in df_motif_list.iterrows():
    df_read_type_frac = pd.read_csv(eval_statsDir + '/' + row['file_name'])
    df_read_type_frac2 = df_read_type_frac.loc[df_read_type_frac['chipseq_norm_signalVal'] > 0].copy()
    df_read_type_frac_all = pd.concat([df_read_type_frac_all, df_read_type_frac2], ignore_index=True)
df_read_type_frac_all['chipseq_norm_signalVal_log10'] = np.log10(df_read_type_frac_all['chipseq_norm_signalVal'])
df_read_type_frac_all_long = df_read_type_frac_all.melt(id_vars=['cell_type', 'TF_name', 'motif_str', 'read_type', 'chipseq_norm_signalVal_log10'], value_vars=['read_frac_defined', 'read_frac_predicted', 'read_frac_predicted_prob'], var_name='read_frac_type', value_name='read_frac_value')
df_read_type_frac_all_long['read_frac_type'] = df_read_type_frac_all_long['read_frac_type'].str.replace('read_frac_', '')

df_spearman_corr = df_read_type_frac_all_long.groupby(['cell_type', 'TF_name', 'motif_str', 'read_frac_type', 'read_type']).apply(lambda grp: pd.Series(spearmanr_wrapper(grp['read_frac_value'], grp['chipseq_norm_signalVal_log10'])), include_groups=False).reset_index()
# df_spearman_corr['pval'] = df_spearman_corr['pval'].fillna(1)  # Replace NaN p-values with 1
# df_spearman_corr['pval_corrected'] = multipletests(df_spearman_corr['pval'], method='fdr_bh')[1]  # p value correction

if summary_eval_statsDir is not None:
    df_auprc_all.to_csv(summary_eval_statsDir + '/auprc_all.csv', index=False)
    df_auroc_all.to_csv(summary_eval_statsDir + '/auroc_all.csv', index=False)
    df_read_type_frac_all.to_csv(summary_eval_statsDir + '/read_type_frac_all.csv', index=False)
    df_spearman_corr.to_csv(summary_eval_statsDir + '/spearman_corr.csv', index=False)

if plot_outDir is not None:
    plt = (
        ggplot(df_f1_all, aes(x='model_type', y='f1_weighted', fill='model_type')) + 
        geom_jitter(width=0.2, stroke=0, alpha=0.5) +
        geom_boxplot(fill=None, width=0.4, outlier_alpha=0.25) +
        theme_light()
    )
    plt.save(plot_outDir + '/f1_weighted_boxplot.pdf', width=5, height=3)

    plt = (
        ggplot(df_auprc_all, aes(x='model_type', y='AUPRC', fill='model_type')) + 
        geom_jitter(width=0.2, stroke=0, alpha=0.5) +
        geom_boxplot(fill=None, width=0.4, outlier_alpha=0.25) + 
        facet_wrap('read_type') + theme_light()
    )
    plt.save(plot_outDir + '/auprc_boxplot.pdf', width=11, height=3)

    plt = (
        ggplot(df_auroc_all, aes(x='model_type', y='AUROC', fill='model_type')) + 
        geom_jitter(width=0.2, stroke=0, alpha=0.5) +
        geom_boxplot(fill=None, width=0.4, outlier_alpha=0.25) + 
        facet_wrap('read_type') + theme_light()
    )
    plt.save(plot_outDir + '/auroc_boxplot.pdf', width=11, height=3)

    plt = (
        ggplot(df_spearman_corr, aes(x='read_frac_type', y='rho', fill='read_frac_type')) + 
        geom_jitter(width=0.2, stroke=0, alpha=0.5) +
        geom_boxplot(fill=None, width=0.4, outlier_alpha=0.25) +
        facet_wrap('~read_type') +
        ggtitle('Spearman_corr: read_frac ~ chipseq_norm_signalVal_log10') +
        theme_light()
    )
    plt.save(plot_outDir + '/spearman_corr_boxplot.pdf', width=12, height=4)

# filter motif stats file
if (motif_stats_path is not None) and (filtered_motif_stats_path is not None) and (cell_type is not None):
    df_motif_stats = pd.read_csv(motif_stats_path)
    df_spearman_corr2 = df_spearman_corr.loc[(df_spearman_corr['read_type']=='bound') & (df_spearman_corr['read_frac_type'] == 'predicted_prob')].copy()
    selected_TFs = df_spearman_corr2.loc[(df_spearman_corr2['rho'].abs()>0.2) & (df_spearman_corr2['pval']<0.05), 'TF_name'].unique()
    df_motif_stats = df_motif_stats.loc[(df_motif_stats['TF'].isin(selected_TFs)) & (df_motif_stats['cell_type']==cell_type)].copy()
    # df_motif_stats = df_motif_stats.loc[(df_motif_stats['cell_type']==cell_type)].copy()

    df_temp = df_spearman_corr2[['cell_type', 'TF_name', 'rho', 'pval']].copy()
    df_temp.rename(columns={'TF_name': 'TF', 'rho': 'model_rho', 'pval': 'model_pval'}, inplace=True)
    df_motif_stats = pd.merge(df_motif_stats, df_temp, on=['TF', 'cell_type'], how='left')
    df_motif_stats.to_csv(filtered_motif_stats_path, index=False)

