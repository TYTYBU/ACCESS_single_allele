import os
import argparse
import pyranges as pr
import pandas as pd
import numpy as np
import pysam
import gzip
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script preprocess motif bed files for downstream ACCESS analysis. It outputs a bed containing the active (bound) and inactive (unbound) motifs with the read counts. ',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Required parameters
    parser.add_argument('factorbookDir', 
                        type=str, 
                        default=None,
                        help=('Factorbook motifs directory. \n'
                              'Default: None'))

    parser.add_argument('k562_bam_in', 
                        type=str, 
                        default=None,
                        help=('K562 BAM file containing the aligned reads. \n'
                              'Default: None'))

    parser.add_argument('hepg2_bam_in', 
                        type=str, 
                        default=None,
                        help=('HepG2 BAM file containing the aligned reads. \n'
                              'Default: None'))

    parser.add_argument('outDir',
                        type=str,
                        default=None,
                        help=('Bed file output directory. \n'
                              'Default: None'))

    parser.add_argument('--start_idx',
                        type=int,
                        default=0,
                        help=('Index in factorbook_MEME_motif_summary_with_chipseq_peaks.csv to start processing motifs. \n'
                              'Default: 0'))

    parser.add_argument('--end_idx',
                        type=int,
                        default=-1,
                        help=('Index in factorbook_MEME_motif_summary_with_chipseq_peaks.csv to end processing motifs. \n'
                              'Default: -1 (the last motif)'))

    parser.add_argument('--flank_len',
                        type=int,
                        default=150,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 150'))

    parser.add_argument('--read_ct_thres',
                        type=int,
                        default=50,
                        help=('Minimum number of reads need to align to a motif, \n'
                              'Default: 50'))

    parser.add_argument('--stats_out',
                        type=str,
                        default=None,
                        help=('Output motif statistics in csv format.'))

    return parser.parse_args()


def get_flank_pos(df, flank_len=300):
    df['flank_start'] = df['start'] - flank_len
    df['flank_end'] = df['end'] + flank_len
    df['center'] = df.apply(lambda row: int((row['start']+row['end'])/2) if row['strand']=='+' else int((row['start']+row['end'])/2+0.5), axis=1)
    return(df)

def motif_row_2_read_ct(bam, row, cover_center=False):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']
    TF_center = row['center']

    read_ct = 0
    for read in bam.fetch(chr, flank_start, flank_end):
        if cover_center:
            if read.reference_start <= TF_center < read.reference_end:
                read_ct += 1
        else:        
            read_ct += 1   
    return(read_ct)


# load arguments
args = parse_args()
factorbook_dir = args.factorbookDir
bam_path_k562 = args.k562_bam_in
bam_path_hepg2 = args.hepg2_bam_in
our_dir = args.outDir
start_idx = args.start_idx
end_idx = args.end_idx
flank_len = args.flank_len
read_ct_thres = args.read_ct_thres
stats_out_path = args.stats_out


# load bam file
logging.info(f'Loading K562 bam file from {bam_path_k562}')
if not os.path.exists(bam_path_k562 + '.bai'):
    logging.info(f'Creating index for bam file.')
    pysam.index(bam_path_k562)

logging.info(f'Loading HepG2 bam file from {bam_path_hepg2}')
if not os.path.exists(bam_path_hepg2 + '.bai'):
    logging.info(f'Creating index for bam file.')
    pysam.index(bam_path_hepg2)

bamList = {'K562': pysam.AlignmentFile(bam_path_k562,'rb'), 
           'HepG2': pysam.AlignmentFile(bam_path_hepg2,'rb')}


# process motifs from factorbook
logging.info(f'Loading motifs from {factorbook_dir}.')
df = pd.read_csv(factorbook_dir + 'factorbook_MEME_motif_summary_with_chipseq_peaks.csv', index_col=0)

if start_idx < 0:
    start_idx = 0
elif start_idx >= len(df):
    start_idx = len(df)-1

if end_idx < 0:
    end_idx = len(df)-1
elif end_idx >= len(df):
    end_idx = len(df)-1
logging.info(f'Processing motifs from index {start_idx} to {end_idx}.')

os.makedirs(our_dir + '/K562', exist_ok=True)
os.makedirs(our_dir + '/HepG2', exist_ok=True)
logging.info(f'Output motif files to {our_dir}.')

tfs, cell_types, motifs, active_motif_ct, inactive_motif_ct = [], [], [], [], []
for i in range(start_idx, end_idx+1):
    tf, cell_type, identifier, motif = df['Target'][i], df['Biosample'][i], df['Identifier'][i], df['motif'][i]

    # create grs for chipseq peaks
    chip_peak_grs = pr.read_bed(f'{factorbook_dir}/chipseq_peaks/{cell_type}/{tf}_{identifier}.bed.gz')

    # read motif matching results from fimo
    motif_match_file = f'{factorbook_dir}/motif_matching/{cell_type}/{tf}_{motif}/fimo.tsv'
    df_motif = pd.read_csv(motif_match_file, sep='\t', on_bad_lines='skip')

    if(len(df_motif) >= 3):  
        df_motif = df_motif.drop(columns=['motif_alt_id'])
        df_motif = df_motif.dropna(axis=0)
        df_motif['start'] = df_motif['start'].astype(int)
        df_motif['stop'] = df_motif['stop'].astype(int)
        df_motif = df_motif[df_motif['sequence_name'].isin(chip_peak_grs.Chromosome.unique())]
        df_motif['fimo_pVal'] = -np.log10(df_motif['p-value'])

        # create grs for mached motifs
        motif_match_grs = pr.from_dict({'Chromosome': df_motif['sequence_name'],
                                        'Start': df_motif['start'],
                                        'End': df_motif['stop'],
                                        'Name': tf,
                                        'Score': df_motif['score'],
                                        'Strand': df_motif['strand'],
                                        'fimo_pVal': df_motif['fimo_pVal']})
        motif_match_grs = motif_match_grs.sort()

        # overlap with ChIP-seq peaks to get labels
        motif_match_tp_grs = motif_match_grs.overlap(chip_peak_grs, strandedness=False, invert=False)
        motif_match_tn_grs = motif_match_grs.overlap(chip_peak_grs, strandedness=False, invert=True)

        # add chipseq score to motifs
        motif_match_tp_grs_added_score = motif_match_tp_grs.join(chip_peak_grs, how="left").drop(like='_b').drop(['ThickEnd', 'BlockCount'])
        df_active = motif_match_tp_grs_added_score.df
        df_active = df_active.rename(columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end', 'Name': 'TF',
                                        'Score': 'score', 'Strand': 'strand', 'ThickStart': 'chipseq_signalVal', 'ItemRGB': 'chipseq_qVal'})
        
        motif_match_tn_grs_added_score = motif_match_tn_grs.join(chip_peak_grs, how="left").drop(like='_b').drop(['ThickEnd', 'BlockCount'])
        df_inactive = motif_match_tn_grs_added_score.df
        df_inactive = df_inactive.rename(columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end', 'Name': 'TF',
                                        'Score': 'score', 'Strand': 'strand', 'ThickStart': 'chipseq_signalVal', 'ItemRGB': 'chipseq_qVal'})
        df_active['motif_type'] = 'Active'
        df_inactive['motif_type'] = 'Inactive'

        # add read count, filter for motifs with at least 50 reads
        bam = bamList[cell_type]

        df_active = get_flank_pos(df_active, flank_len)
        df_active['motif_read_ct'] = df_active.apply(lambda row: motif_row_2_read_ct(bam, row, cover_center=True), axis=1)
        df_active2 = df_active.loc[df_active['motif_read_ct']>=read_ct_thres].sort_values('motif_read_ct', ascending=False).copy()

        df_inactive = get_flank_pos(df_inactive, flank_len)
        df_inactive['motif_read_ct'] = df_inactive.apply(lambda row: motif_row_2_read_ct(bam, row, cover_center=True), axis=1)
        df_inactive2 = df_inactive.loc[df_inactive['motif_read_ct']>=read_ct_thres].sort_values('motif_read_ct', ascending=False).copy()

        df_motif_out2 = pd.concat([df_active2, df_inactive2])
        df_motif_out2.to_csv(f'{our_dir}/{cell_type}/{tf}_{motif}.motifs.csv.gz', index=False, compression='gzip')
        
        tfs.append(tf)
        motifs.append(motif)
        cell_types.append(cell_type)
        active_motif_ct.append(len(df_active2))
        inactive_motif_ct.append(len(df_inactive2))


if args.stats_out is not None:
    df_stats = pd.DataFrame(data={'cell_type': cell_types, 'TF': tfs, 'motif': motifs, 
                                  'active_motif_ct': active_motif_ct, 'inactive_motif_ct': inactive_motif_ct})
    df_stats.to_csv(stats_out_path, index=False)















