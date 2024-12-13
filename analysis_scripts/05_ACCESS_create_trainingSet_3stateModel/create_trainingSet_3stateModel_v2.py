import os
import argparse
import pandas as pd
import numpy as np
import pysam
import logging
from access_util import *
from plotnine import *
import h5py
import ast

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def parse_args():
    parser = argparse.ArgumentParser(
        description='This script calculates TF motif features and then use the TF specific features to label reads from bound (active) motifs into bound, unbound and recently bound. ',
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
    
    parser.add_argument('h5Dir', 
                        type=str, 
                        default=None,
                        help=('Output directory for the training dataset (h5ad format). \n'
                              'Default: None'))
    
    parser.add_argument('--predefined_labels_outDir',
                        type=str,
                        default=None,
                        help=('Output directory for edit fraction data.\n'
                              'Default: None'))

    parser.add_argument('--flank_len',
                        type=int,
                        default=150,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 150'))
    
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
featuresDir = args.featuresDir

cell_type = args.cell_type
TF_name = args.TF_name
motif_str = args.motif_str

h5Dir = args.h5Dir
predefined_labels_outDir = args.predefined_labels_outDir

flank_len = args.flank_len
max_read_ct_label = args.max_read_ct_label


# load feature file
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
bed_path = featuresDir + f'/clean_motifs/{cell_type}_{TF_name}_{motif_str}.clean_motifs.csv.gz'
logging.info(f'Loading clean motif bed file from {bed_path}')
df_clean_motifs = pd.read_csv(bed_path)
df_active_ref = df_clean_motifs.loc[(df_clean_motifs['motif_type'] == 'Active') & (df_clean_motifs['use_type'] == 'train')].copy()
df_active_noRef = df_clean_motifs.loc[(df_clean_motifs['motif_type'] == 'Active') & (df_clean_motifs['use_type'] == 'test')].copy()
logging.info(f'Active motifs for training: {len(df_active_ref)}, for testing: {len(df_active_noRef)}')


# create training dataset
logging.info(f'Creating training set ...')
oneHotMat_dict = df_motif_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, df_active_ref,  
                                                         center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                                         trough_thres, peak_thres, global_thres,
                                                         flank_len = flank_len, read_ct_thres = max_read_ct_label)
x_train = np.vstack((oneHotMat_dict['unbound'], oneHotMat_dict['bound'], oneHotMat_dict['recently bound']))
y_train_class = [0] * oneHotMat_dict['unbound'].shape[0] + [1] * oneHotMat_dict['bound'].shape[0] + [2] * oneHotMat_dict['recently bound'].shape[0]

# calculate class weights
samples_per_class = np.array([oneHotMat_dict['unbound'].shape[0], oneHotMat_dict['bound'].shape[0], oneHotMat_dict['recently bound'].shape[0]])
class_weights = sum(samples_per_class) / (len(samples_per_class) * samples_per_class)
logging.info(f'Class weights: {class_weights}')

# create test dataset
logging.info(f'Creating testing set ...')
oneHotMat_dict = df_motif_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, df_active_noRef,  
                                                         center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                                         trough_thres, peak_thres, global_thres,
                                                         flank_len = flank_len, read_ct_thres = max_read_ct_label)
x_test = np.vstack((oneHotMat_dict['unbound'], oneHotMat_dict['bound'], oneHotMat_dict['recently bound']))
y_test_class = [0] * oneHotMat_dict['unbound'].shape[0] + [1] * oneHotMat_dict['bound'].shape[0] + [2] * oneHotMat_dict['recently bound'].shape[0]

# save training data to disk
if h5Dir is not None:
    h5_path = h5Dir + f'/{cell_type}_{TF_name}_{motif_str}__3readTypes_oneHot.h5'
    logging.info(f'Writing training set to: {h5_path}')
    hf = h5py.File(h5_path, 'w')
    hf.create_dataset('x_train', data=x_train, compression="gzip")
    hf.create_dataset('y_train_class', data=y_train_class, compression="gzip")
    hf.create_dataset('class_weights', data=class_weights, compression="gzip")
    hf.create_dataset('x_test', data=x_test, compression="gzip")
    hf.create_dataset('y_test_class', data=y_test_class, compression="gzip")
    hf.close()

    df_clean_motifs.to_hdf(h5_path, key='df_clean_motifs', mode='a')


if predefined_labels_outDir is not None:
    if not os.path.exists(f'{predefined_labels_outDir}/plots/'):
        os.makedirs(f'{predefined_labels_outDir}/plots/')

    x_train = np.array(x_train)
    y_train_class = np.array(y_train_class)

    logging.info(f'calculating edit fractions for pre-defined read labels ...')
    df_editFrac_train = oneHotMat_2_edit_frac_by_read_type_labels(x_train, y_train_class)
    df_trainSet_stats = pd.DataFrame({'Unbound': [oneHotMat_dict['unbound'].shape[0]], 'Bound': [oneHotMat_dict['bound'].shape[0]], 'Recently bound': [oneHotMat_dict['recently bound'].shape[0]]})
    
    logging.info(f'Output pre-defined read labels data to {predefined_labels_outDir}')
    df_editFrac_train.to_csv(predefined_labels_outDir + f'/{cell_type}_{TF_name}_{motif_str}.predefined_labels_editFrac.csv', index=False)
    df_trainSet_stats.to_csv(predefined_labels_outDir + f'/{cell_type}_{TF_name}_{motif_str}.predefined_labels_ct.csv', index=False)
    
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

    figure_path = predefined_labels_outDir + f'/plots/{cell_type}_{TF_name}_{motif_str}.predefined_labels.pdf'
    logging.info(f'Output edit fraction plot for pre-defeined labelled reads to: {figure_path}')
    plt.save(filename=figure_path, width=7, height=4)
















