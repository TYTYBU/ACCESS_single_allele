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
    parser.add_argument('stats_summary_in_path', 
                        type=str, 
                        default=None,
                        help=('Selected motif summary stats file input path. \n'
                              'Default: None'))

    parser.add_argument('featuresDir', 
                        type=str, 
                        default=None,
                        help=('Motif feature part files directory. \n'
                              'Default: None'))

    parser.add_argument('stats_outDir', 
                        type=str, 
                        default=None,
                        help=('Motif summary stats file with feature data output directiory. \n'
                              'Default: None'))
    
    parser.add_argument('--spearman_thres', 
                        type=float, 
                        default=0.1,
                        help=('Minimum spearman correlation rho for motif pass filtering. \n'
                              'Default: 0.1'))

    return parser.parse_args()


# load arguments
args = parse_args()
stats_summary_in_path = args.stats_summary_in_path
featuresDir = args.featuresDir
stats_outDir = args.stats_outDir
spearman_thres = args.spearman_thres


logging.info(f'Loading motif stats summary from: {stats_summary_in_path}')
df_stats_all = pd.read_csv(stats_summary_in_path)

df_motif_features_all = pd.DataFrame()
for index, row in df_stats_all.iterrows():
    cell_type, TF_name, motif_str = row['cell_type'], row['TF'], row['motif']
    
    features_path = featuresDir + f'/{cell_type}_{TF_name}_{motif_str}.motif_features.csv'
    if os.path.exists(features_path):
        df_motif_features = pd.read_csv(features_path)
    else:
        logging.warning(f'Feature file not found: {features_path}')
        continue

    df_spearman_corr = pd.read_csv(f'/net/bgm/sherwood/factorbook_data/UltimaGen/motif_features_v2/weights_selection/spearman_corr/{cell_type}_{TF_name}_{motif_str}.spearman_corr_all_weights.csv')
    df_spearman_corr_bound = df_spearman_corr.loc[(df_spearman_corr['read_type'] == 'bound') & (df_spearman_corr['pval']<0.1)].copy()
    best_rho = df_spearman_corr_bound['rho'].max()
    best_weights = df_spearman_corr_bound.loc[df_spearman_corr_bound['rho'] == best_rho].agg({'peak_weight': 'mean', 'footprint_weight': 'mean', 'global_std_weight': 'mean', 'rho': 'mean', 'pval': 'mean'})

    df_motif_features['peak_weight'] = best_weights['peak_weight']
    df_motif_features['footprint_weight'] = best_weights['footprint_weight']
    df_motif_features['global_std_weight'] = best_weights['global_std_weight']

    df_motif_features[['cell_type', 'TF', 'motif']] = [cell_type, TF_name, motif_str]
    df_motif_features_all = pd.concat([df_motif_features_all, df_motif_features])

df_stats_all = pd.merge(df_stats_all, df_motif_features_all, on=['cell_type', 'TF', 'motif'], how='left')
df_stats_all = df_stats_all.dropna(subset=['rho'])

df_stats_filtered = df_stats_all.loc[(df_stats_all['rho'] >= spearman_thres) & (df_stats_all['pval'] < 0.05)].copy()
df_stats_filtered_HepG2 = df_stats_filtered.loc[df_stats_filtered['cell_type'] == 'HepG2'].copy()
df_stats_filtered_K562 = df_stats_filtered.loc[df_stats_filtered['cell_type'] == 'K562'].copy()

if stats_outDir is not None:
    logging.info(f'Output motif stats summary with feature info to: {stats_outDir}')
    df_stats_all.to_csv(stats_outDir + '/TF_summary.features_all.csv', index=False)
    df_stats_filtered.to_csv(stats_outDir + '/TF_summary.features_filtered.csv', index=False)
    df_stats_filtered_HepG2.to_csv(stats_outDir + '/TF_summary.features_filtered_HepG2.csv', index=False)
    df_stats_filtered_K562.to_csv(stats_outDir + '/TF_summary.features_filtered_K562.csv', index=False)











