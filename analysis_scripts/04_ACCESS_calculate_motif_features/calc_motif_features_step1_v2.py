import os
import argparse
import pandas as pd
import numpy as np
import pysam
import gzip
import logging
from access_util import *
from plotnine import *

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='This script calculates TF motif features. ',
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

    parser.add_argument('bedDir', 
                        type=str, 
                        default=None,
                        help=('Motif bed files directory. \n'
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

    parser.add_argument('features_outDir', 
                        type=str, 
                        default=None,
                        help=('Output TF motif feature values in csv format. \n'
                              'Default: None'))

    parser.add_argument('--fimo_pVal_thres',
                        type=float,
                        default=0.0001,
                        help=('Min fimo score for both active and inactive motifs. \n'
                              'Default: 0.0001'))

    parser.add_argument('--chipseq_qVal_thres',
                        type=float,
                        default=0.05,
                        help=('Min chipseq qVal for active motifs. For inactive motifs: -1. \n'
                              'Default: 0.05'))

    parser.add_argument('--flank_len',
                        type=int,
                        default=150,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 150'))

    parser.add_argument('--sampled_motif_ct',
                        type=int,
                        default=200,
                        help=('Number of motifs sampled to calculate composite edit fractions. \n'
                              'Default: 200'))
    
    parser.add_argument('--max_read_ct_label',
                        type=int,
                        default=1000000,
                        help=('Max number of reads to store in memory and process. \n'
                              'Default: 1M'))

    return parser.parse_args()


# load arguments
args = parse_args()
bam_path = args.bamFile
fasta_path = args.fastaFile
bedDir = args.bedDir
cell_type = args.cell_type
TF_name = args.TF_name
motif_str = args.motif_str
features_outDir = args.features_outDir
fimo_pVal_thres = -np.log10(args.fimo_pVal_thres)
chipseq_qVal_thres = -np.log10(args.chipseq_qVal_thres)
flank_len = args.flank_len
sampled_motif_ct = args.sampled_motif_ct
max_read_ct_label = args.max_read_ct_label

feature_chr_list = ['chr1']
train_chr_list = ['chr'+str(y) for y in [x for x in range(3,12)] + [x for x in range(13,22)] + ['X', 'Y']]
test_chr_list = ['chr'+str(y) for y in [x for x in [2,12,22]]]


# load bam file
logging.info(f'Loading bam file from {bam_path}')
if not os.path.exists(bam_path + '.bai'):
    logging.info(f'Creating index for bam file.')
    pysam.index(bam_path)
bam = pysam.AlignmentFile(bam_path,'rb')

# load genome file
logging.info(f'Loading genome from {fasta_path}')
if not os.path.exists(fasta_path + '.fai'):
    print("fasta file need to indexed first!")
    exit(1)
fasta = pysam.FastaFile(fasta_path)

# load motif bed file
bed_path = bedDir + f'/{cell_type}/{TF_name}_{motif_str}.motifs.csv.gz'
logging.info(f'Loading motif bed file from {bed_path}')
df_motifs = pd.read_csv(bed_path)

# filter for standard active and inactive motifs
df_motifs_filtered = filter_df_motifs(df_motifs, fimo_pVal_thres=fimo_pVal_thres, chipseq_qVal_thres=chipseq_qVal_thres)
df_motifs_filtered = get_flank_pos(df_motifs_filtered, flank_len)
df_active_pool = df_motifs_filtered.loc[df_motifs_filtered['motif_type'] == 'Active'].copy()
df_inactive_pool = df_motifs_filtered.loc[df_motifs_filtered['motif_type'] == 'Inactive'].copy()
logging.info(f'filtered active motifs: {len(df_active_pool)}, inactive motifs: {len(df_inactive_pool)}')


# select motifs for feature extraction
df_active_features_pool = df_active_pool.loc[df_active_pool['chr'].isin(feature_chr_list)].copy()
df_active_features_pool.sort_values(by=['chipseq_signalVal'], ascending=False, inplace=True)
feature_motif_ct = int(len(df_active_features_pool) * 0.25)
feature_motif_ct = max(min(100, len(df_active_features_pool)), feature_motif_ct)
feature_motif_ct = min(sampled_motif_ct, feature_motif_ct)
df_active_features = df_active_features_pool.head(feature_motif_ct).copy()
df_active_unused = df_active_features_pool.tail(len(df_active_features_pool) - feature_motif_ct).copy()
df_active_features['use_type'] = 'feature'
df_inactive_features_pool = df_inactive_pool.loc[df_inactive_pool['chr'].isin(feature_chr_list)].copy()
df_inactive_features = df_inactive_features_pool.sample(n=feature_motif_ct).copy()
logging.info(f'Active motifs reserved for feature extraction: {len(df_active_features_pool)}, used: {len(df_active_features)}')

# select active motifs for training and testing
df_active_ref = df_active_pool.loc[df_active_pool['chr'].isin(train_chr_list)].copy()
df_active_ref = pd.concat([df_active_ref, df_active_unused], ignore_index=True)
df_active_ref['use_type'] = 'train'
df_active_noRef = df_active_pool.loc[df_active_pool['chr'].isin(test_chr_list)].copy()
df_active_noRef['use_type'] = 'test'
logging.info(f'Active motifs for training: {len(df_active_ref)}, for testing: {len(df_active_noRef)}')

# select inactive motifs for training and testing
df_inactive_ref = df_inactive_pool.loc[df_inactive_pool['chr'].isin(train_chr_list)].copy()
df_inactive_ref['use_type'] = 'train'
df_inactive_noRef = df_inactive_pool.loc[df_inactive_pool['chr'].isin(test_chr_list)].copy()
df_inactive_noRef['use_type'] = 'test'
logging.info(f'Inactive motifs for training: {len(df_inactive_ref)}, for testing: {len(df_inactive_noRef)}')

df_clean_motifs = pd.concat([df_active_features_pool, df_active_ref, df_active_noRef, df_inactive_ref, df_inactive_noRef])


logging.info(f'Calculating motif features ...')
df_editFrac_features = df_motifs_2_df_reads_2_composite_edit_frac(bam, df_active_features, flank_len)
df_editFrac_features['edit_frac_ma'] = df_editFrac_features['edit_frac'].rolling(window=15, min_periods=1, center=True).mean()
df_editFrac_inactive_sampled = df_motifs_2_df_reads_2_composite_edit_frac(bam, df_inactive_features, flank_len)
df_editFrac_inactive_sampled['edit_frac_ma'] = df_editFrac_inactive_sampled['edit_frac'].rolling(window=15, min_periods=1, center=True).mean()
df_motif_features = calc_motif_features_v2(df_editFrac_features, df_editFrac_inactive_sampled, flank_len)
df_active_motif_features_expanded = expand_active_motif_features(df_motif_features)

logging.info(f'Writing clean motif bed file ...')
os.makedirs(features_outDir + '/clean_motifs/', exist_ok=True)
df_clean_motifs.to_csv(features_outDir + f'/clean_motifs/{cell_type}_{TF_name}_{motif_str}.clean_motifs.csv.gz', index=False, compression='gzip')

logging.info(f'Writing bound and unbound motif edit fractions ...')
os.makedirs(features_outDir + '/unlabeled_editFrac/plots/', exist_ok=True)
df_editFrac_features.to_csv(features_outDir + f'/unlabeled_editFrac/{cell_type}_{TF_name}_{motif_str}.edit_frac_features.csv', index=False)
df_editFrac_inactive_sampled.to_csv(features_outDir + f'/unlabeled_editFrac/{cell_type}_{TF_name}_{motif_str}.edit_frac_inactive_sampled.csv', index=False)

plt = (
    ggplot() + 
    geom_rect(df_active_motif_features_expanded, aes(xmin='relative_pos-0.5', xmax='relative_pos+0.5', fill='feature'), ymin=0, ymax=1, alpha=0.25) +
    geom_segment(df_active_motif_features_expanded, aes(x='relative_pos-0.5', xend='relative_pos+0.5', y='mean_editFrac', yend='mean_editFrac', color='feature'), size=1) + 
    geom_line(df_editFrac_features, aes(x='relative_pos', y='edit_frac_ma'), color='red', size=0.5) + 
    # geom_line(df_editFrac_features, aes(x='relative_pos', y='edit_frac'), color='red', size=0.25) + 
    stat_smooth(df_editFrac_features, aes(x='relative_pos', y='edit_frac'), color='red', method='mavg', method_args={'window': 3, 'center': True}, size=0.25, se=False) + 
    geom_line(df_editFrac_inactive_sampled, aes(x='relative_pos', y='edit_frac_ma'), color='blue', size=0.5) + 
    # geom_line(df_editFrac_inactive_sampled, aes(x='relative_pos', y='edit_frac'), color='blue', size=0.25) + 
    stat_smooth(df_editFrac_inactive_sampled, aes(x='relative_pos', y='edit_frac'), color='blue', method='mavg', method_args={'window': 3, 'center': True}, size=0.25, se=False) + 
    scale_color_manual(values=['#e78ac3', '#a6d854', '#ffd92f']) + scale_fill_manual(values=['#e78ac3', '#a6d854', '#ffd92f']) +
    labs(x='Relative position', y='Edit fraction') + ggtitle(f'{cell_type} {TF_name} {motif_str}') +
    theme_light()
)
plt.save(filename=features_outDir + f'/unlabeled_editFrac/plots/{cell_type}_{TF_name}_{motif_str}.edit_frac_features.pdf', width=7, height=4)


# sampling active motifs to label reads
df_active_ref_filtered = df_active_ref.loc[df_active_ref['chipseq_signalVal'] > 0].copy()
df_active_ref_sampled = df_active_ref_filtered.sample(sampled_motif_ct if len(df_active_ref) > sampled_motif_ct else len(df_active_ref))
df_active_ref_sampled['chipseq_norm_signalVal'] = df_active_ref_sampled['chipseq_signalVal'] / df_active_ref_sampled['chipseq_signalVal'].max()

logging.info(f'Iterating through weights combos ...')
df_read_type_frac_params_combo_gene = pd.DataFrame()
for index, row in df_active_ref_sampled.iterrows():
    df_read_type_frac_params_combo_row = get_read_type_frac_for_param_combos(bam, row, df_motif_features, flank_len=flank_len)
    df_read_type_frac_params_combo_row['chipseq_signalVal'] = row['chipseq_signalVal']
    df_read_type_frac_params_combo_row['chipseq_norm_signalVal'] = row['chipseq_norm_signalVal']
    df_read_type_frac_params_combo_gene = pd.concat([df_read_type_frac_params_combo_gene, df_read_type_frac_params_combo_row], ignore_index=True)

# spearman correlation 
read_type_list = ['unbound', 'bound', 'recently bound']
logging.info(f'Calculating spearman correlation for each weights combo ...')
df_read_type_frac_params_combo_gene['read_type'] = pd.Categorical(df_read_type_frac_params_combo_gene['read_type'], categories=read_type_list)
df_read_type_frac_params_combo_gene['chipseq_norm_signalVal_log10'] = np.log10(df_read_type_frac_params_combo_gene['chipseq_norm_signalVal'])
df_read_type_frac_params_combo_gene_long = df_read_type_frac_params_combo_gene.melt(id_vars=['peak_weight', 'footprint_weight', 'global_std_weight', 'read_type', 'read_ct_defined', 'chipseq_norm_signalVal_log10'], value_vars='read_frac_defined', var_name='read_frac_type', value_name='read_frac_value')
df_read_type_frac_params_combo_gene_long['read_frac_type'] = df_read_type_frac_params_combo_gene_long['read_frac_type'].str.replace('read_frac_', '')
df_spearman_corr = df_read_type_frac_params_combo_gene_long.groupby(['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'], observed=False).apply(lambda grp: pd.Series(spearmanr_wrapper(grp['read_frac_value'], grp['chipseq_norm_signalVal_log10'])), include_groups=False).reset_index()

df_read_ct = df_read_type_frac_params_combo_gene_long.groupby(['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'], observed=False).agg(read_ct=('read_ct_defined', 'sum')).reset_index()
df_spearman_corr = pd.merge(df_spearman_corr, df_read_ct, on=['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'])
df_spearman_corr = df_spearman_corr.sort_values(by=['read_type', 'rho'], ascending=[True, False])
df_spearman_corr_bound = df_spearman_corr.loc[(df_spearman_corr['read_type'] == 'bound') & (df_spearman_corr['pval']<0.1)].copy()

logging.info(f'Writing read type fraction data and correlation data ...')
os.makedirs(features_outDir + '/weights_selection/read_type_frac', exist_ok=True)
df_read_type_frac_params_combo_gene.to_csv(features_outDir + f'/weights_selection/read_type_frac/{cell_type}_{TF_name}_{motif_str}.read_type_frac_all_weights.csv.gz', index=False, compression='gzip')

os.makedirs(features_outDir + '/weights_selection/spearman_corr/plots', exist_ok=True)
df_spearman_corr.to_csv(features_outDir + f'/weights_selection/spearman_corr/{cell_type}_{TF_name}_{motif_str}.spearman_corr_all_weights.csv', index=False)

if len(df_spearman_corr_bound) > 0:
    # generate correlation plots for read type fractions
    best_rho = df_spearman_corr_bound['rho'].max()
    best_weights = df_spearman_corr_bound.loc[df_spearman_corr_bound['rho'] == best_rho].agg({'peak_weight': 'mean', 'footprint_weight': 'mean', 'global_std_weight': 'mean', 'rho': 'mean', 'pval': 'mean'})
    df_read_type_frac_params_combo_gene_long_best = df_read_type_frac_params_combo_gene_long.loc[(df_read_type_frac_params_combo_gene_long['peak_weight'] == best_weights['peak_weight']) & (df_read_type_frac_params_combo_gene_long['footprint_weight'] == best_weights['footprint_weight']) & (df_read_type_frac_params_combo_gene_long['global_std_weight'] == best_weights['global_std_weight'])].copy()
    df_spearman_corr_best = df_spearman_corr.loc[(df_spearman_corr['peak_weight'] == best_weights['peak_weight']) & (df_spearman_corr['footprint_weight'] == best_weights['footprint_weight']) & (df_spearman_corr['global_std_weight'] == best_weights['global_std_weight'])].copy()
    df_spearman_corr_best['spearman_str'] = df_spearman_corr_best.apply(lambda row: f'spearman rho={row['rho']:.2f}\np={row['pval']:.2e}', axis=1)
    df_temp = df_read_type_frac_params_combo_gene_long_best.groupby(['read_frac_type', 'read_type']).apply(lambda grp: pd.Series({'x': grp['read_frac_value'].max(), 'y': grp['chipseq_norm_signalVal_log10'].max()}), include_groups=False).reset_index()
    df_spearman_corr_best = pd.merge(df_spearman_corr_best, df_temp, on=['read_frac_type', 'read_type'])

    plt = (
        ggplot() + 
        geom_point(df_read_type_frac_params_combo_gene_long_best, aes(x='read_frac_value', y='chipseq_norm_signalVal_log10', color='chipseq_norm_signalVal_log10'), alpha=0.25, size=0.25) + 
        scale_color_gradient(low="blue", high="red") + 
        geom_text(df_spearman_corr_best, aes(label='spearman_str', x='x', y='y'), ha='right', va='top') + 
        theme_light() + facet_grid('read_frac_type ~ read_type', scales='free') + 
        labs(x='Read fractions in motifs', y='-log10 normalized ChIP-Seq signal', color='-log10 normalized ChIP-Seq signal') + 
        ggtitle(f'{cell_type} {TF_name} {motif_str}') +
        theme(legend_position='bottom')
    )
    plt.save(filename=features_outDir + f'/weights_selection/spearman_corr/plots/{cell_type}_{TF_name}_{motif_str}.spearman_corr_best.pdf', width=8, height=4)

    logging.info(f'Writing motif feature stats ...')
    df_motif_features = pd.concat([df_motif_features, pd.DataFrame([best_weights[['peak_weight', 'footprint_weight', 'global_std_weight', 'rho', 'pval']]]).reset_index(drop=True)], axis=1)
    df_motif_features[df_spearman_corr_best['read_type']] = df_spearman_corr_best['read_ct']
    df_motif_features.to_csv(f'{features_outDir}/{cell_type}_{TF_name}_{motif_str}.motif_features.csv', index=False)


    # labeled read edit fraction plot from best weights
    center_trough_pos_list = df_motif_features['center_pos_list'][0]
    left_peak_pos_list = df_motif_features['l_peak_pos_list'][0]
    right_peak_pos_list = df_motif_features['r_peak_pos_list'][0]

    peak_w = df_motif_features['peak_weight'][0]
    footprint_w = df_motif_features['footprint_weight'][0]
    global_std_w = df_motif_features['global_std_weight'][0]
    logging.info(f'Best weights combo:')
    logging.info(f'Footprint weight: {footprint_w}, Peak weight: {peak_w}, Global std weight: {global_std_w}')

    peak_thres = df_motif_features['left_peak'][0]*peak_w*0.5 + df_motif_features['right_peak'][0]*peak_w*0.5 + df_motif_features['footprint'][0]*(1-peak_w)
    trough_thres = df_motif_features['left_peak'][0]*footprint_w*0.5 + df_motif_features['right_peak'][0]*footprint_w*0.5 + df_motif_features['footprint'][0]*(1-footprint_w)
    global_thres = df_motif_features['inactive_mean'][0] + df_motif_features['inactive_std'][0]*global_std_w

    logging.info(f'Creating sample oneHotMat dictionary ...')
    oneHotMat_dict = df_motif_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, df_active_ref_sampled,  
                                                             center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                                             trough_thres, peak_thres, global_thres,
                                                             flank_len = flank_len, read_ct_thres = max_read_ct_label)
    x_train = np.vstack((oneHotMat_dict['unbound'], oneHotMat_dict['bound'], oneHotMat_dict['recently bound']))
    y_train_class = [0] * oneHotMat_dict['unbound'].shape[0] + [1] * oneHotMat_dict['bound'].shape[0] + [2] * oneHotMat_dict['recently bound'].shape[0]
    x_train = np.array(x_train)
    y_train_class = np.array(y_train_class)

    logging.info(f'calculating edit fractions for pre-defined read labels ...')
    df_editFrac_train = oneHotMat_2_edit_frac_by_read_type_labels(x_train, y_train_class)
    df_trainSet_stats = pd.DataFrame({'Unbound': [oneHotMat_dict['unbound'].shape[0]], 'Bound': [oneHotMat_dict['bound'].shape[0]], 'Recently bound': [oneHotMat_dict['recently bound'].shape[0]]})
    
    logging.info(f'Output labeled edit fraction data ...')
    os.makedirs(features_outDir + '/labeled_editFrac/plots/', exist_ok=True)
    df_editFrac_train.to_csv(features_outDir + f'/labeled_editFrac/{cell_type}_{TF_name}_{motif_str}.predefined_labels_editFrac.csv', index=False)
    df_trainSet_stats.to_csv(features_outDir + f'/labeled_editFrac/{cell_type}_{TF_name}_{motif_str}.predefined_labels_ct.csv', index=False)
    
    df_trainSet_stats['ct_str'] = df_trainSet_stats.apply(lambda row: f"Unbound: {row['Unbound']}\nBound: {row['Bound']}\nRecently bound: {row['Recently bound']}", axis=1)
    df_trainSet_stats['motif_str'] = motif_str
    df_active_motif_features_expanded['motif_str'] = motif_str
    df_editFrac_train['motif_str'] = motif_str

    plt = (
        ggplot() + 
        geom_label(df_trainSet_stats, aes(label='ct_str'), y=0.8, x=-150, color='black', fill='white', size=8, ha='left') +
        geom_rect(df_active_motif_features_expanded, aes(xmin='relative_pos-0.5', xmax='relative_pos+0.5', fill='feature'), ymin=0, ymax=1, alpha=0.25) +
        geom_line(df_editFrac_train, aes(x='relative_pos', y='edit_frac', color='read_type'), size=0.2) +
        stat_smooth(df_editFrac_train, aes(x='relative_pos', y='edit_frac', color='read_type'), method='mavg', method_args={'window': 3, 'center': True}, size=0.5, se=False) + 
        scale_color_manual(values=['red', 'green', 'blue']) + scale_fill_manual(values=['#e78ac3', '#a6d854', '#ffd92f']) +
        labs(x='Relative position', y='Edit fraction') + facet_wrap('motif_str') + ylim(0, 1) + 
        ggtitle(f"{cell_type} {TF_name} pre-defined label edit fraction") + 
        theme_light()
    )

    figure_path = features_outDir + f'/labeled_editFrac/plots/{cell_type}_{TF_name}_{motif_str}.predefined_labels.pdf'
    logging.info(f'Output edit fraction plot for pre-defeined labeled reads to: {figure_path}')
    plt.save(filename=figure_path, width=7, height=4)
else:
    logging.info(f'No weights combo achieved significant spearman correlation')
    best_weights = pd.Series({'peak_weight': np.nan, 'footprint_weight': np.nan, 'global_std_weight': np.nan, 
                              'rho': np.nan, 'pval': np.nan,
                              'unbound': np.nan, 'bound': np.nan, 'recently bound': np.nan})
    df_motif_features = pd.concat([df_motif_features, pd.DataFrame([best_weights]).reset_index(drop=True)], axis=1)
    df_motif_features.to_csv(f'{features_outDir}/{cell_type}_{TF_name}_{motif_str}.motif_features.csv', index=False)




# df_read_type_frac_params_combo_gene = pd.read_csv(f'{features_outDir}/weights_selection/read_type_frac/{cell_type}_{TF_name}_{motif_str}.read_type_frac_all_weights.csv.gz')
# read_type_list = ['unbound', 'bound', 'recently bound']

# df_read_type_frac_params_combo_gene['read_type'] = pd.Categorical(df_read_type_frac_params_combo_gene['read_type'], categories=read_type_list)
# df_read_type_frac_params_combo_gene['chipseq_norm_signalVal_log10'] = np.log10(df_read_type_frac_params_combo_gene['chipseq_norm_signalVal'])
# df_read_type_frac_params_combo_gene_long = df_read_type_frac_params_combo_gene.melt(id_vars=['peak_weight', 'footprint_weight', 'global_std_weight', 'read_type', 'read_ct_defined', 'chipseq_norm_signalVal_log10'], value_vars='read_frac_defined', var_name='read_frac_type', value_name='read_frac_value')
# df_read_type_frac_params_combo_gene_long['read_frac_type'] = df_read_type_frac_params_combo_gene_long['read_frac_type'].str.replace('read_frac_', '')
# df_spearman_corr = df_read_type_frac_params_combo_gene_long.groupby(['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'], observed=False).apply(lambda grp: pd.Series(spearmanr_wrapper(grp['read_frac_value'], grp['chipseq_norm_signalVal_log10'])), include_groups=False).reset_index()
# df_read_ct = df_read_type_frac_params_combo_gene_long.groupby(['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'], observed=False).agg(read_ct=('read_ct_defined', 'sum')).reset_index()
# df_spearman_corr = pd.merge(df_spearman_corr, df_read_ct, on=['peak_weight', 'footprint_weight', 'global_std_weight', 'read_frac_type', 'read_type'])
# df_spearman_corr = df_spearman_corr.sort_values(by=['read_type', 'rho'], ascending=[True, False])

# best_weights = df_spearman_corr.loc[df_spearman_corr['read_type'] == 'unbound'].iloc[0]
# df_read_type_frac_params_combo_gene_long_best = df_read_type_frac_params_combo_gene_long.loc[(df_read_type_frac_params_combo_gene_long['peak_weight'] == best_weights['peak_weight']) & (df_read_type_frac_params_combo_gene_long['footprint_weight'] == best_weights['footprint_weight']) & (df_read_type_frac_params_combo_gene_long['global_std_weight'] == best_weights['global_std_weight'])].copy()
# df_spearman_corr_best = df_spearman_corr.loc[(df_spearman_corr['peak_weight'] == best_weights['peak_weight']) & (df_spearman_corr['footprint_weight'] == best_weights['footprint_weight']) & (df_spearman_corr['global_std_weight'] == best_weights['global_std_weight'])].copy()
# df_spearman_corr_best['spearman_str'] = df_spearman_corr_best.apply(lambda row: f'spearman rho={row['rho']:.2f}\np={row['pval']:.2e}', axis=1)
# df_temp = df_read_type_frac_params_combo_gene_long_best.groupby(['read_frac_type', 'read_type']).apply(lambda grp: pd.Series({'x': grp['read_frac_value'].max(), 'y': grp['chipseq_norm_signalVal_log10'].max()}), include_groups=False).reset_index()
# df_spearman_corr_best = pd.merge(df_spearman_corr_best, df_temp, on=['read_frac_type', 'read_type'])
# df_motif_features = pd.read_csv(f'{features_outDir}/{cell_type}_{TF_name}_{motif_str}.motif_features.csv')
# df_motif_features = df_motif_features[['center_pos_list', 'l_peak_pos_list', 'r_peak_pos_list', 'footprint', 'left_peak', 'right_peak']].copy()
# df_motif_features = pd.concat([df_motif_features, pd.DataFrame([best_weights[['peak_weight', 'footprint_weight', 'global_std_weight', 'rho', 'pval']]]).reset_index(drop=True)], axis=1)
# df_motif_features[df_spearman_corr_best['read_type']] = df_spearman_corr_best['read_ct']
# df_motif_features.to_csv(f'{features_outDir}/{cell_type}_{TF_name}_{motif_str}.motif_features.csv', index=False)













