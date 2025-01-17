import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from plotnine import *

from scipy.fft import fft
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

from access_util import *


rho_thres = 0.2
read_ct_thres = 100

one_side_window_size = 1
proximity_range = [16, 100]
bound_filter_str_list = ['_filtered_001_099'] 
motif_pair_orientation_list = ['either strand']
adj_method = ['mean of footprint boundaries', 
              'mean of left peak + footprint + right peak boundaries',
              'median of footprint positions',
              'mean of footprint positions',
              'median of left peak + footprint + right peak positions',
              'mean of left peak + footprint + right peak positions']
adj_method_num = 3
exp_colname = 'expected_probability'
oe_colname = 'delta_oe_prob'
oe_colname_lite = 'delta_oe'
period_range = [8, 13]

df_cobinding_bound_by_pos_all = pd.DataFrame()
df_cobinding_bound_ma_all = pd.DataFrame()
df_fft_all = pd.DataFrame()
df_fft_stats_all = pd.DataFrame()
for cell_type in ['HepG2', 'K562']:
    # motif stats file
    statsFile = f'./motif_stats/TF_summary.features_filtered_{cell_type}_modeled.csv'
    df_motif_stats = pd.read_csv(statsFile)
    df_motif_stats['center_pos_list'] = df_motif_stats['center_pos_list'].apply(lambda x: ast.literal_eval(x))
    df_motif_stats['l_peak_pos_list'] = df_motif_stats['l_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
    df_motif_stats['r_peak_pos_list'] = df_motif_stats['r_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
    df_motif_stats = df_motif_stats.loc[(df_motif_stats['cell_type'] == cell_type) & (df_motif_stats['model_rho'] > rho_thres)].copy()

    # cobinding stats file
    cobind_statsDir = f'./occuPIE_co_occupancy/{cell_type}/'
    stats_file_suffix = '.cobinding_stats.csv.gz'
    df_motif_list = pd.DataFrame([x.replace(stats_file_suffix, '').split('_', maxsplit=2) + [x] for x in os.listdir(cobind_statsDir) if stats_file_suffix in x], columns=['cell_type', 'TF_name', 'motif_str', 'file_name'])

    # iterate through each TF
    for TF1_name in tqdm(df_motif_list['TF_name'].unique(), desc=f"Processing {cell_type} TFs"):
        # load new data
        row_cobind = df_motif_list.loc[df_motif_list['TF_name'] == TF1_name].iloc[0].copy()
        df_cobinding = pd.read_csv(cobind_statsDir + '/' + row_cobind['file_name'])
        df_cobinding_bound = df_cobinding.loc[
            (df_cobinding['read_type_pair'] == 'bound--bound') & 
            (df_cobinding['shared_read_ct'] >= read_ct_thres) &
            (df_cobinding['TF2'].isin(df_motif_list['TF_name']))
        ].copy()

        df_cobinding_bound['delta_oe_prob'] = df_cobinding_bound['observed_probability'] - df_cobinding_bound['expected_probability']

        for bound_filter_str in bound_filter_str_list:
            if bound_filter_str == '':
                bound_filter = False
            elif bound_filter_str == '_filtered_001_099':
                bound_filter = True
                bound_filter_range = [0.01, 0.99]
            elif bound_filter_str == '_filtered_010_090':
                bound_filter = True
                bound_filter_range = [0.10, 0.90]

            # filtering out motif pairs with extremely low or high bound probability on either motifs
            if bound_filter:
                df_cobinding_bound_filtered = df_cobinding_bound.loc[df_cobinding_bound[exp_colname].between(bound_filter_range[0], bound_filter_range[1], inclusive='both')].copy()
            else:
                df_cobinding_bound_filtered = df_cobinding_bound

            # adjust center distance calculation
            df_cobinding_bound_filtered['center_dist'] = df_cobinding_bound_filtered.apply(lambda row: calc_center_distance_v2(row, df_motif_stats, adj_method_num=adj_method_num), axis=1)

            for motif_pair_orientation in motif_pair_orientation_list:
                # differentiate TF1 TF2 motif pair orientation
                if motif_pair_orientation == 'same strand':
                    df_cobinding_bound_temp = df_cobinding_bound_filtered.loc[df_cobinding_bound_filtered['strand1'] == df_cobinding_bound_filtered['strand2']].copy()
                elif motif_pair_orientation == 'opposite strand':
                    df_cobinding_bound_temp = df_cobinding_bound_filtered.loc[df_cobinding_bound_filtered['strand1'] != df_cobinding_bound_filtered['strand2']].copy()
                else:
                    df_cobinding_bound_temp = df_cobinding_bound_filtered

                if len(df_cobinding_bound_temp) > 0:
                    # normalize observed - expected values
                    filtered_values = df_cobinding_bound_temp[oe_colname].replace([np.inf, -np.inf], np.nan).dropna()
                    oe_perc05 = np.percentile(filtered_values, 5)
                    oe_perc95 = np.percentile(filtered_values, 95)
                    df_cobinding_bound_temp[oe_colname + '_norm'] = df_cobinding_bound_temp[oe_colname].apply(lambda x: (x - oe_perc05)/(oe_perc95 - oe_perc05))

                    # calculate median observed - expected values for each position
                    df_cobinding_bound_by_pos = agg_by_bin_window(df_cobinding_bound_temp, bin_col=oe_colname + '_norm', max_center_dist=proximity_range[1], one_side_window_size=0)
                    df_cobinding_bound_by_pos = df_cobinding_bound_by_pos.loc[df_cobinding_bound_by_pos['center_dist'].abs() >= proximity_range[0]].copy()
                    df_cobinding_bound_by_pos['cell_type'] = cell_type
                    df_cobinding_bound_by_pos['TF_name'] = TF1_name
                    df_cobinding_bound_by_pos['bound_filter'] = str(bound_filter_range) if bound_filter else 'no bound filter'
                    df_cobinding_bound_by_pos['motif_pair_orientation'] = motif_pair_orientation
                    df_cobinding_bound_by_pos.rename(columns={oe_colname + '_norm' + '_by_pos': oe_colname_lite + '_by_pos' + '_norm'}, inplace=True)
                    df_cobinding_bound_by_pos_all = pd.concat([df_cobinding_bound_by_pos_all, df_cobinding_bound_by_pos], ignore_index=True)

                    # FFT analysis with randomization test
                    df_fft = df_cobinding_2_df_fft_v2(df_cobinding_bound_temp, y_column=oe_colname + '_norm', y_column_rename=oe_colname_lite + '_by_pos' + '_norm', 
                                                    proximity_range = [16, 100], one_side_window_size=0, 
                                                    randomization_test=True, num_iterations=2000)
                    
                    df_fft['cell_type'] = cell_type
                    df_fft['TF_name'] = TF1_name
                    df_fft['bound_filter'] = str(bound_filter_range) if bound_filter else 'no bound filter'
                    df_fft['motif_pair_orientation'] = motif_pair_orientation
                    df_fft_all = pd.concat([df_fft_all, df_fft], ignore_index=True)


outDir = f'./peiodicity_strength/'
os.makedirs(outDir, exist_ok=True)

# p-value correction
df_fft_filtered_all = df_fft_all.loc[
    (df_fft_all['bound_filter'] == '[0.01, 0.99]') &
    (df_fft_all['motif_pair_orientation'] == 'either strand')
].copy()

df_fft_stats_all = pd.DataFrame()
for (cell_type, TF_name), df_fft in df_fft_filtered_all.groupby(['cell_type', 'TF_name']):
    df_fft_stats = get_dominant_period(df_fft, min_period=period_range[0], max_period=period_range[1], alpha=1)
    df_fft_stats['pval_bh'] = multipletests(df_fft_stats['pval'], method='fdr_bh')[1]
    df_fft_stats['pval_bonferroni'] = multipletests(df_fft_stats['pval'], method='bonferroni')[1]

    df_fft_stats['cell_type'] = cell_type
    df_fft_stats['TF_name'] = TF_name
    df_fft_stats_all = pd.concat([df_fft_stats_all, df_fft_stats], ignore_index=True)


df_fft_all.to_csv(f'{outDir}/TF1-all.fft.csv.gz', index=False, compression='gzip')
df_cobinding_bound_by_pos_all.to_csv(f'{outDir}/TF1-all.cobinding_bound_by_pos.csv.gz', index=False, compression='gzip')
df_fft_stats_all.to_csv(f'{outDir}/TF1-all.fft_stats_all.csv', index=False)
                            