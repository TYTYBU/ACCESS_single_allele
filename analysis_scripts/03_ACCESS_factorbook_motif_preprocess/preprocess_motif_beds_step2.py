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
        description='This script summarizes and filters motifs based on various stats. For TFs with multiple motif bed files, the representitive motif is selected by e-value and then by active (bound) motif count. ',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('factorbook_stats_in', 
                        type=str, 
                        default=None,
                        help=('Motif stats file from factorbook in csv. It contains e-value for each motif. \n'
                              'Default: None'))
    
    parser.add_argument('bedDir', 
                        type=str, 
                        default=None,
                        help=('Motif bed files directory. \n'
                              'Default: None'))

    parser.add_argument('stats_in_prefix', 
                        type=str, 
                        default=None,
                        help=('Motif stats input file prefix. \n'
                              'Default: None'))

    parser.add_argument('stats_summary_out_prefix',
                        type=str,
                        default=None,
                        help=('Output directory for combined motif statistics in csv format.'))

    parser.add_argument('start_pt',
                        type=int,
                        default=None,
                        help=('Staring part number for motif stats files. \n'
                              'Default: None'))

    parser.add_argument('end_pt',
                        type=int,
                        default=None,
                        help=('Ending part number for motif stats files. \n'
                              'Default: None'))
    
    parser.add_argument('--motif_list_in',
                        type=str,
                        default=None,
                        help=('A list of TF motif name (one per line) to include in the filtering process. \n'
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

    parser.add_argument('--active_motif_ct_thres',
                        type=int,
                        default=1000,
                        help=('Minimum active motif count threshold. \n'
                              'Default: 1000'))

    return parser.parse_args()

def select_motif(df_grp):
    # Sort by e_value (ascending) and filtered_active_motif_ct (descending)
    return df_grp.sort_values(by=['e_value', 'filtered_active_motif_ct'], ascending=[True, False]).head(1)

# load arguments
args = parse_args()
factorbook_stats_in = args.factorbook_stats_in
bedDir = args.bedDir
stats_in_prefix = args.stats_in_prefix
stats_summary_out_prefix = args.stats_summary_out_prefix
start_pt = args.start_pt
end_pt = args.end_pt
motif_list_in = args.motif_list_in
fimo_pVal_thres = -np.log10(args.fimo_pVal_thres)
chipseq_qVal_thres = -np.log10(args.chipseq_qVal_thres)
active_motif_ct_thres = args.active_motif_ct_thres


logging.info(f'Reading TF motif bed files from {bedDir}')
df_stats_all = pd.DataFrame()
for pt in range(start_pt, end_pt+1):
    stats_path = stats_in_prefix + f'.pt{pt}.csv'
    df_stats = pd.read_csv(stats_path)

    filtered_active_motif_ct = []
    filtered_inactive_motif_ct = []
    for index, row in df_stats.iterrows():
        cell_type, TF_name, motif_str = row['cell_type'], row['TF'], row['motif']
        bed_path = bedDir + f'/{cell_type}/{TF_name}_{motif_str}.motifs.csv.gz'
        df_motifs = pd.read_csv(bed_path)

        df_active_pool = df_motifs.loc[(df_motifs['fimo_pVal'] >= fimo_pVal_thres) & ((df_motifs['motif_type'] == 'Active') & (df_motifs['chipseq_qVal'] > chipseq_qVal_thres))].copy()
        df_inactive_pool = df_motifs.loc[(df_motifs['fimo_pVal'] >= fimo_pVal_thres) & ((df_motifs['motif_type'] == 'Inactive') & (df_motifs['chipseq_qVal'] < 0))].copy()

        # df_active_pool = df_motifs.loc[(df_motifs['score'] >= 15) & ((df_motifs['motif_type'] == 'Active') & (df_motifs['chipseq_signalVal'] >= 50))].copy()
        # df_inactive_pool = df_motifs.loc[(df_motifs['score'] >= 15) & ((df_motifs['motif_type'] == 'Inactive') & (df_motifs['chipseq_signalVal'] < 0))].copy()

        filtered_active_motif_ct.append(len(df_active_pool))
        filtered_inactive_motif_ct.append(len(df_inactive_pool))

    df_stats['filtered_active_motif_ct'] = filtered_active_motif_ct
    df_stats['filtered_inactive_motif_ct'] = filtered_inactive_motif_ct
    df_stats_all = pd.concat([df_stats_all, df_stats], ignore_index=True)

# select the TF motif in each cell line based on the following criteria:
# 1. included in the master TF list
# 2. have the higher evalue
# 3. have higher active motif ct

logging.info(f'Loading factorbook motif stats from {factorbook_stats_in}')
df_factorbook = pd.read_csv(factorbook_stats_in)
df_factorbook['e_value'] = df_factorbook['e_value'].apply(lambda x: 0 if x == '< 1e-300' else float(x))
df_factorbook.rename(columns={'Target':'TF', 'Biosample':'cell_type'}, inplace=True)
df_factorbook = df_factorbook[['cell_type', 'TF', 'motif', 'e_value']].copy()
df_stats_all = pd.merge(df_stats_all, df_factorbook, on=['cell_type', 'TF', 'motif'], how='left')

df_stats_selected = df_stats_all.groupby(['cell_type', 'TF']).apply(select_motif).reset_index(drop=True)
if factorbook_stats_in is not None:
    logging.info(f'Loading motif list to include from {factorbook_stats_in}')
    df_TF_master = pd.read_csv(motif_list_in)
    df_stats_selected = df_stats_selected.loc[df_stats_selected['TF'].isin(df_TF_master['TF'])].copy()
    
df_stats_selected = df_stats_selected.loc[df_stats_selected['filtered_active_motif_ct'] >= active_motif_ct_thres].reset_index(drop=True)
logging.info(f'TF motifs: {len(df_stats_all)}, passed all filters: {len(df_stats_selected)}')

stats_out_path = stats_summary_out_prefix + '.all.csv'
df_stats_all.to_csv(stats_out_path, index=False)

stats_out_path = stats_summary_out_prefix + '.selected.csv'
df_stats_selected.to_csv(stats_out_path, index=False)














