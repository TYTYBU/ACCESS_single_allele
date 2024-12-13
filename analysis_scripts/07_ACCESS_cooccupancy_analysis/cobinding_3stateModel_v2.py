import os
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from scipy.stats import fisher_exact, mannwhitneyu
from access_util import *
from plotnine import *
import logging
import pickle
import ast

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='Co-binding analysis with defined labels.',
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

    parser.add_argument('statsFile', 
                        type=str, 
                        default=None,
                        help=('Motif feature stats CSV file. \n'
                              'Default: None'))
    
    parser.add_argument('bedDir', 
                        type=str, 
                        default=None,
                        help=('Path to motif bed file directory. \n'
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

    parser.add_argument('outDir',
                        type=str,
                        default=None,
                        help=('Output directory for cobinding stat files. \n'
                              'Default: None'))

    parser.add_argument('--min_dist',
                        type=int,
                        default=5,
                        help=('Minimum distance between motif centers to be considered in proximity. \n'
                              'Default: 5'))
    
    parser.add_argument('--max_dist',
                        type=int,
                        default=100,
                        help=('Maximum distance between motif centers to be considered in proximity. \n'
                              'Default: 100'))

    parser.add_argument('--flank_len',
                        type=int,
                        default=150,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 150'))

    parser.add_argument('--read_ct_thres',
                        type=int,
                        default=50,
                        help=('Minimum number of mapped reads for a motif locus to be eligible for cobinding analysis. \n'
                              'Default: 50'))
    
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
    
    parser.add_argument('--correlation_thres',
                        type=float,
                        default=0.2,
                        help=('Min spearman rho between ChIPseq signal and bound fraction for a motif to be included. \n'
                              'Default: 0.2'))

    parser.add_argument('--make_plots',
                        action='store_true',
                        default=False,
                        help='Generate co-binding plots. \n'
                             'Default: False')

    parser.add_argument('--adj_method_num',
                        type=int,
                        default=3,
                        help=('Methods to adjust motif center positions \n'
                              '0: mean of footprint boundaries\n'
                              '1: mean of left peak + footprint + right peak boundaries\n'
                              '2: mean of footprint positions\n'
                              '3: median of footprint positions\n'
                              '4: mean of left peak + footprint + right peak positions\n'
                              '5: median of left peak + footprint + right peak positions\n'
                              'Default: 3'))
    
    return parser.parse_args()


# load arguments
args = parse_args()
bam_path = args.bamFile
fasta_path=args.fastaFile
statsFile = args.statsFile
bedDir = args.bedDir
modelDir = args.modelDir
cell_type = args.cell_type
TF_name = args.TF_name
motif_str = args.motif_str
outDir = args.outDir
min_dist = args.min_dist
max_dist = args.max_dist
flank_len = args.flank_len
read_ct_thres = args.read_ct_thres
fimo_pVal_thres = -np.log10(args.fimo_pVal_thres)
chipseq_qVal_thres = -np.log10(args.chipseq_qVal_thres)
rho_thres = args.correlation_thres
make_plots = args.make_plots
adj_method_num = args.adj_method_num
adj_method = ['mean of footprint boundaries', 
              'mean of left peak + footprint + right peak boundaries',
              'median of footprint positions',
              'mean of footprint positions',
              'median of left peak + footprint + right peak positions',
              'mean of left peak + footprint + right peak positions']
logging.info(f'Motif center adjustment method: {adj_method[adj_method_num]}')

# parameters
one_side_window_size = 1
search_window = [-x for x in range(1, 9)] + [x for x in range(1, 9)] # how far in adjacent motifs to serach for proximal motifs
proximity_range = [min_dist, max_dist] # proximity range for cobinding motifs
read_type_list = ['unbound', 'bound', 'recently bound']

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

# load all models
logging.info(f'Loading all {cell_type} 3-state model from {modelDir} ...')
model_file_suffix = '__3stateModels.pkl'
df_model_list = pd.DataFrame([x.replace(model_file_suffix, '').split('_', maxsplit=2) + [x] for x in os.listdir(modelDir) if model_file_suffix in x], columns=['cell_type', 'TF_name', 'motif_str', 'file_name'])
models_dict_all = dict()
for index, row in df_model_list.iterrows():
    pkl_path = modelDir + f'/{row['file_name']}'
    with open(pkl_path, 'rb') as file:
        models_dict = pickle.load(file)
        models_dict_all[row['TF_name']] = models_dict

if TF_name not in models_dict_all:
    logging.info(f'Model not found for {TF_name}, co-binding analysis terminated.')
    exit(1)

# load motif stats file
logging.info(f'Loading motifs combined stats file: {statsFile} ...')
df_motif_stats = pd.read_csv(statsFile)
df_motif_stats['center_pos_list'] = df_motif_stats['center_pos_list'].apply(lambda x: ast.literal_eval(x))
df_motif_stats['l_peak_pos_list'] = df_motif_stats['l_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
df_motif_stats['r_peak_pos_list'] = df_motif_stats['r_peak_pos_list'].apply(lambda x: ast.literal_eval(x))
df_motif_stats = df_motif_stats.loc[(df_motif_stats['cell_type'] == cell_type) & (df_motif_stats['model_rho'] > rho_thres)].copy()


# load motif1 bed file
motif1_bed_path = bedDir + f'/{cell_type}/{TF_name}_{motif_str}.motifs.csv.gz'
logging.info(f'Loading motif1 bed file from {motif1_bed_path} ...')
df_motif1 = pd.read_csv(motif1_bed_path).drop_duplicates()
df_motif1 = filter_df_motifs(df_motif1, fimo_pVal_thres=fimo_pVal_thres, chipseq_qVal_thres=chipseq_qVal_thres, read_ct_thres=read_ct_thres)
df_motif1 = get_flank_pos(df_motif1, flank_len)
df_motif1 = df_motif1.sort_values(['TF', 'chr', 'start', 'end', 'strand', 'chipseq_signalVal'], ascending=[True, True, True, True, True, False]).groupby(['TF', 'chr', 'start', 'end', 'strand']).first().reset_index()
df_motif1.loc[df_motif1['chipseq_signalVal']<0, 'chipseq_signalVal'] = 0
df_motif1['chipseq_norm_signalVal'] = df_motif1['chipseq_signalVal'] / df_motif1['chipseq_signalVal'].max()
df_motif1['is_target_motif'] = True

row_TF1 = df_motif_stats.loc[(df_motif_stats['TF'] == TF_name) & (df_motif_stats['motif'] == motif_str)].iloc[0]
df_motif1['center_corrected'] = df_motif1.apply(lambda row: adjust_motif_center_v2(row, row_TF1), axis=1)
df_motif1_active = df_motif1.loc[df_motif1['motif_type'] == 'Active'].copy()


## find active motif1 and active motif2 loci pairs in proximity for motif1 TFs
df_valid_loci_pairs_all = pd.DataFrame()
df_cobinding_all = pd.DataFrame()
logging.info(f'Searching all loci pairs in proximity for {TF_name} motifs ...')
for index, row_TF2 in df_motif_stats.iterrows():
    logging.info(f'Searching for {TF_name} - {row_TF2['TF']} loci pairs  ...')
    motif2_bed_path = bedDir + f'/{cell_type}/{row_TF2['TF']}_{row_TF2['motif']}.motifs.csv.gz'
    df_motif2 = pd.read_csv(motif2_bed_path).drop_duplicates()
    df_motif2 = filter_df_motifs(df_motif2, fimo_pVal_thres=fimo_pVal_thres, chipseq_qVal_thres=chipseq_qVal_thres, read_ct_thres=read_ct_thres)
    df_motif2 = get_flank_pos(df_motif2, flank_len)
    df_motif2 = df_motif2.sort_values(['TF', 'chr', 'start', 'end', 'strand', 'chipseq_signalVal'], ascending=[True, True, True, True, True, False]).groupby(['TF', 'chr', 'start', 'end', 'strand']).first().reset_index()
    df_motif2_active = df_motif2.loc[df_motif2['motif_type'] == 'Active'].copy()
    df_motif2_active.loc[df_motif2_active['chipseq_signalVal']<0, 'chipseq_signalVal'] = 0
    df_motif2_active['chipseq_norm_signalVal'] = df_motif2_active['chipseq_signalVal'] / df_motif2_active['chipseq_signalVal'].max()
    df_motif2_active['is_target_motif'] = False
    df_motif2_active['center_corrected'] = df_motif2_active.apply(lambda row2: adjust_motif_center_v2(row2, row_TF2), axis=1)
    
    df_loci_mixed = pd.concat([df_motif1, df_motif2_active]).sort_values(['chr', 'center_corrected']) # contains loci pairs for all active and inactive motif1
    df_loci_mixed_original = df_loci_mixed.add_suffix('1')
    df_valid_loci_pairs_expanded = pd.DataFrame() # contains all loci pairs < max center distance threshold
    for i in search_window:
        df_loci_mixed_shifted = df_loci_mixed.shift(i, fill_value=0)
        df_loci_mixed_shifted = df_loci_mixed_shifted.add_suffix('2')

        df_loci_pairs = pd.concat([df_loci_mixed_original, df_loci_mixed_shifted], axis=1)
        df_loci_pairs = df_loci_pairs.loc[
            (df_loci_pairs['chr1'] == df_loci_pairs['chr2']) & 
            (df_loci_pairs['is_target_motif1'] == True) & 
            (df_loci_pairs['is_target_motif2'] == False)
        ].copy()

        if len(df_loci_pairs) == 0:
            continue

        # df_loci_pairs['center_dist'] = df_loci_pairs.apply(lambda row: calc_center_distance(row, row_TF1, row_TF2), axis=1)
        df_loci_pairs['center_dist'] = df_loci_pairs.apply(lambda row: calc_center_distance_v2(row, df_motif_stats, adj_method_num=adj_method_num), axis=1)
        # only filter for loci pairs with less than the max center distance threshold
        df_loci_pairs = df_loci_pairs.loc[~ (df_loci_pairs['center_dist'].abs() > proximity_range[1])].copy()
        df_valid_loci_pairs_expanded = pd.concat([df_valid_loci_pairs_expanded, df_loci_pairs], ignore_index=True)

    df_valid_loci_pairs_expanded.drop_duplicates(inplace=True, ignore_index=True)
    df_valid_loci_pairs_expanded_acitveOnly = df_valid_loci_pairs_expanded.loc[df_valid_loci_pairs_expanded['motif_type1'] == 'Active'].copy()
    df_valid_loci_pairs = df_valid_loci_pairs_expanded.loc[df_valid_loci_pairs_expanded['center_dist'].abs().between(proximity_range[0]-one_side_window_size, proximity_range[1], inclusive='both')].copy() # further filter for loci pairs with more than the min center distance threshold
    df_valid_loci_pairs_acitveOnly = df_valid_loci_pairs.loc[df_valid_loci_pairs['motif_type1'] == 'Active'].copy()
    df_valid_loci_pairs_all = pd.concat([df_valid_loci_pairs_all, df_valid_loci_pairs_acitveOnly], ignore_index=True)

    if row_TF2['TF'] not in models_dict_all:
        logging.info(f'Model not found for {row_TF2['TF']}. Skipping ...')
        continue
    
    # calculate observed - expected for loci pairs
    logging.info(f'Calculating expected - observed read for motifs in proximity ...')
    df_cobinding_motif2 = df_loci_pair_2_df_cobinding_predicted(bam, fasta, df_valid_loci_pairs_acitveOnly, df_motif_stats, models_dict_all, flank_len=flank_len, shared_read_ct_thres=read_ct_thres)
    df_cobinding_motif2['read_type_pair'] = df_cobinding_motif2['read_type_1'] + '--' + df_cobinding_motif2['read_type_2'] 
    df_cobinding_motif2['delta_oe_frac'] = df_cobinding_motif2['observed_read_fraction'] - df_cobinding_motif2['expected_read_fraction']
    df_cobinding_motif2['delta_oe_prob'] = df_cobinding_motif2['observed_probability'] - df_cobinding_motif2['expected_probability']
    df_cobinding_all = pd.concat([df_cobinding_all, df_cobinding_motif2], ignore_index=True)
    
    if make_plots:
        # generate allele level position-resolved plots for individual motif2
        os.makedirs(outDir + f'/plots/TF1_TF2/{cell_type}_{TF_name}/', exist_ok=True)
        df_cobinding_motif2_bound = df_cobinding_motif2.loc[df_cobinding_motif2['read_type_pair'] == 'bound--bound'].copy()

        df_cobinding_motif2_bound_by_pos = agg_by_bin_window(df_cobinding_motif2_bound, bin_col='delta_oe_prob', max_center_dist=max_dist, one_side_window_size=0)
        df_cobinding_motif2_bound_by_pos = df_cobinding_motif2_bound_by_pos.loc[df_cobinding_motif2_bound_by_pos['center_dist'].abs() >= proximity_range[0]].copy()
        df_cobinding_motif2_bound_by_pos['read_type_pair'] = 'bound--bound'
        df_cobinding_motif2_bound_ma = agg_by_bin_window(df_cobinding_motif2_bound, bin_col='delta_oe_prob', max_center_dist=max_dist, one_side_window_size=one_side_window_size)
        df_cobinding_motif2_bound_ma = df_cobinding_motif2_bound_ma.loc[df_cobinding_motif2_bound_ma['center_dist'].abs() >= proximity_range[0]].copy()
        df_cobinding_motif2_bound_ma['read_type_pair'] = 'bound--bound'

        plt = (
            ggplot() + 
            geom_point(df_cobinding_motif2_bound_by_pos, aes(x='center_dist', y='delta_oe_prob_by_pos'), color='blue') +
            geom_line(df_cobinding_motif2_bound_ma, aes(x='center_dist', y='delta_oe_prob_by_pos'), color='purple', size=0.5) +
            labs(x='Distance to adjacent motif', y='Observed - Expected probability') + 
            facet_wrap('read_type_pair') + ggtitle(f'{TF_name} - {row_TF2['TF']} in proximity') + theme_light() + 
            theme(figure_size=(12,4), legend_position='bottom')
        )
        vline_positions = list(range(-99, 100, 11))
        for pos in vline_positions:
            plt += geom_vline(xintercept=pos, linetype='dashed', color='red', size=0.25)
        out_path = outDir + f'/plots/TF1_TF2/{cell_type}_{TF_name}/{TF_name}_{row_TF2['TF']}.delta_oe_prob_by_pos.pdf'
        plt.save(out_path)
        logging.info(f'Output TF1-TF2 positional observed-expected probability plots to {out_path}')

        df_cobinding_motif2_bound_by_pos = agg_by_bin_window(df_cobinding_motif2_bound, bin_col='delta_oe_frac', max_center_dist=max_dist, one_side_window_size=0)
        df_cobinding_motif2_bound_by_pos = df_cobinding_motif2_bound_by_pos.loc[df_cobinding_motif2_bound_by_pos['center_dist'].abs() >= proximity_range[0]].copy()
        df_cobinding_motif2_bound_by_pos['read_type_pair'] = 'bound--bound'
        df_cobinding_motif2_bound_ma = agg_by_bin_window(df_cobinding_motif2_bound, bin_col='delta_oe_frac', max_center_dist=max_dist, one_side_window_size=one_side_window_size)
        df_cobinding_motif2_bound_ma = df_cobinding_motif2_bound_ma.loc[df_cobinding_motif2_bound_ma['center_dist'].abs() >= proximity_range[0]].copy()
        df_cobinding_motif2_bound_ma['read_type_pair'] = 'bound--bound'

        plt = (
            ggplot() + 
            geom_point(df_cobinding_motif2_bound_by_pos, aes(x='center_dist', y='delta_oe_frac_by_pos'), color='blue') +
            geom_line(df_cobinding_motif2_bound_ma, aes(x='center_dist', y='delta_oe_frac_by_pos'), color='purple', size=0.5) +
            labs(x='Distance to adjacent motif', y='Observed - Expected fraction') + 
            facet_wrap('read_type_pair') + ggtitle(f'{TF_name} - {row_TF2['TF']} in proximity') + theme_light() + 
            theme(figure_size=(12,4), legend_position='bottom')
        )
        vline_positions = list(range(-99, 100, 11))
        for pos in vline_positions:
            plt += geom_vline(xintercept=pos, linetype='dashed', color='red', size=0.25)
        out_path = outDir + f'/plots/TF1_TF2/{cell_type}_{TF_name}/{TF_name}_{row_TF2['TF']}.delta_oe_frac_by_pos.pdf'
        plt.save(out_path)
        logging.info(f'Output TF1-TF2 positional observed-expected fraction plots to {out_path}')


logging.info(f'Saving all output files to {outDir} ...')
df_valid_loci_pairs_all.to_csv(outDir + f'/{cell_type}_{TF_name}_{motif_str}.loci_pairs.csv.gz', index=False, compression='gzip')
df_cobinding_all.to_csv(outDir + f'/{cell_type}_{TF_name}_{motif_str}.cobinding_stats.csv.gz', index=False, compression='gzip')

if make_plots:
    # generate allel level position-resolved plot for all motif2 combined
    os.makedirs(outDir + f'/plots/TF1_all/', exist_ok=True)
    df_cobinding_bound = df_cobinding_all.loc[df_cobinding_all['read_type_pair'] == 'bound--bound'].copy()

    df_cobinding_bound_by_pos = agg_by_bin_window(df_cobinding_bound, bin_col='delta_oe_prob', max_center_dist=max_dist, one_side_window_size=0)
    df_cobinding_bound_by_pos = df_cobinding_bound_by_pos.loc[df_cobinding_bound_by_pos['center_dist'].abs() >= proximity_range[0]].copy()
    df_cobinding_bound_by_pos['read_type_pair'] = 'bound--bound'
    df_cobinding_bound_ma = agg_by_bin_window(df_cobinding_bound, bin_col='delta_oe_prob', max_center_dist=max_dist, one_side_window_size=one_side_window_size)
    df_cobinding_bound_ma = df_cobinding_bound_ma.loc[df_cobinding_bound_ma['center_dist'].abs() >= proximity_range[0]].copy()
    df_cobinding_bound_ma['read_type_pair'] = 'bound--bound'

    plt = (
        ggplot() + 
        geom_point(df_cobinding_bound_by_pos, aes(x='center_dist', y='delta_oe_prob_by_pos'), color='blue') +
        geom_line(df_cobinding_bound_ma, aes(x='center_dist', y='delta_oe_prob_by_pos'), color='purple', size=0.5) +
        labs(x='Distance to adjacent motif', y='Observed - Expected probability') + 
        facet_wrap('read_type_pair') + ggtitle(f'{TF_name} - all motifs in proximity') + theme_light() + 
        theme(figure_size=(12,4), legend_position='bottom')
    )
    vline_positions = list(range(-99, 100, 11))
    for pos in vline_positions:
        plt += geom_vline(xintercept=pos, linetype='dashed', color='red', size=0.25)
    out_path = outDir + f'/plots/TF1_all/{cell_type}_{TF_name}_{motif_str}.delta_oe_prob_by_pos.pdf'
    plt.save(out_path)
    logging.info(f'TF1-all positional observed-expected probability plots to: {out_path}')

    df_cobinding_bound_by_pos = agg_by_bin_window(df_cobinding_bound, bin_col='delta_oe_frac', max_center_dist=max_dist, one_side_window_size=0)
    df_cobinding_bound_by_pos = df_cobinding_bound_by_pos.loc[df_cobinding_bound_by_pos['center_dist'].abs() >= proximity_range[0]].copy()
    df_cobinding_bound_by_pos['read_type_pair'] = 'bound--bound'
    df_cobinding_bound_ma = agg_by_bin_window(df_cobinding_bound, bin_col='delta_oe_frac', max_center_dist=max_dist, one_side_window_size=one_side_window_size)
    df_cobinding_bound_ma = df_cobinding_bound_ma.loc[df_cobinding_bound_ma['center_dist'].abs() >= proximity_range[0]].copy()
    df_cobinding_bound_ma['read_type_pair'] = 'bound--bound'

    plt = (
        ggplot() + 
        geom_point(df_cobinding_bound_by_pos, aes(x='center_dist', y='delta_oe_frac_by_pos'), color='blue') +
        geom_line(df_cobinding_bound_ma, aes(x='center_dist', y='delta_oe_frac_by_pos'), color='purple', size=0.5) +
        labs(x='Distance to adjacent motif', y='Observed - Expected fraction') + 
        facet_wrap('read_type_pair') + ggtitle(f'{TF_name} - all motifs in proximity') + theme_light() + 
        theme(figure_size=(12,4), legend_position='bottom')
    )
    vline_positions = list(range(-99, 100, 11))
    for pos in vline_positions:
        plt += geom_vline(xintercept=pos, linetype='dashed', color='red', size=0.25)
    out_path = outDir + f'/plots/TF1_all/{cell_type}_{TF_name}_{motif_str}.delta_oe_frac_by_pos.pdf'
    plt.save(out_path)
    logging.info(f'TF1-all positional observed-expected fraction plots to: {out_path}')

