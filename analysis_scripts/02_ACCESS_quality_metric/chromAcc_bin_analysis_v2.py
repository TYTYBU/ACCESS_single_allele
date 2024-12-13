import os
import argparse
import pandas as pd
import numpy as np
import pysam
import pyBigWig
from Bio.Seq import Seq
from plotnine import *
from plotnineseqsuite.logo import geom_logo
from plotnineseqsuite.theme import theme_seq
import logging

import warnings
warnings.filterwarnings("ignore")

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script calculates edit fractions for different chromatin accessibility bins.',
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
                        help=('Fasta file of the genome. \n'
                              'Default: None'))

    parser.add_argument('bwFile', 
                        type=str, 
                        default=None,
                        help=('Accessibility bigwig file from ENCODE. \n'
                              'Default: None'))

    parser.add_argument('--acc_bins', 
                        type=str, 
                        default='0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25',
                        help=('Accessibility bins separated by comma. \n'
                              'Default: 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25'))
        
    parser.add_argument('--baseCt_grouped',
                        type=str,
                        default=None,
                        help=('Output composite edit fraction for each motif and cell line as CSV. \n'
                              'Default: None'))
    
    parser.add_argument('--no_ref_n',
                        action="store_true", 
                        default=False,
                        help=('Do not separate trinucleotide edit motifs. \n'
                              'Default: False'))

    parser.add_argument('--segment_len',
                        type=int,
                        default=1000000,
                        help=('Number of processing segments in each chromosome. \n'
                              'Default: 1000000'))

    parser.add_argument('--explore_mode',
                        action="store_true", 
                        default=False,
                        help=('Counting all possible base changes. Do not calculate edit fractions. Do not convert orientations for trinucleotide edit motifs. \n'
                              'Default: False'))

    return parser.parse_args()


# load arguments
args = parse_args()
bam_path=args.bamFile
fasta_path=args.fastaFile
bw_path=args.bwFile
segment_len = args.segment_len
bin_str = args.acc_bins

if args.explore_mode:
    logging.info(f'Explore mode on. Counting all possible base changes.')

# load bam file
logging.info(f'Loading bam file from {bam_path}')
if not os.path.exists(bam_path + '.bai'):
    logging.info(f'Creating index for bam file.')
    pysam.index(bam_path)
bam = pysam.AlignmentFile(bam_path,'rb')

# load genome file
logging.info(f'Loading genome from {fasta_path}')
fasta = pysam.FastaFile(fasta_path)

# load accessibility score bigwig file
logging.info(f'Loading bigwig file from {bw_path}')
bw = pyBigWig.open(bw_path)


# chromain accessibility bin analysis
df_baseCt_grouped = pd.DataFrame()
motif_span = 1 if not args.no_ref_n else 0
valid_chrom = set(fasta.references) & set(bw.chroms().keys())

for i in range(fasta.nreferences):
    chrom = fasta.references[i]

    if chrom in valid_chrom:
        max_pos = fasta.lengths[i]
        segments = int((max_pos-motif_span*2)/segment_len)+1 # leave motif_span space for the very begnning and the very end of the 
        logging.info(f'{chrom} total number of segments: {segments}')
    else:
        logging.info(f'skipping {chrom} ...')
        continue

    df_acc = pd.DataFrame(bw.intervals(chrom), columns=['start', 'end', 'score'])
    bins = [0] + [float(x) for x in bin_str.split(',')] + [float('inf')]
    df_acc['bin'] = pd.cut(df_acc['score'], bins=bins, right=False)
    df_acc['group'] = ((df_acc['start'] > df_acc['end'].shift()) | (df_acc['bin'] != df_acc['bin'].shift())).cumsum()
    df_acc=df_acc.groupby('group').agg({'start':'min', 'end':'max', 'bin':'first'}).reset_index(drop=True)

    for s in range(0, segments):
        segment_start_pos = motif_span+segment_len*s
        segment_end_pos = min(segment_start_pos+segment_len, max_pos-motif_span)
        logging.info(f'Processing {chrom}, segment {s}, start_pos: {segment_start_pos}, end_pos: {segment_end_pos}')
        ref_seq = fasta.fetch(chrom, segment_start_pos, segment_end_pos)

        (A_ct, C_ct, G_ct, T_ct) = bam.count_coverage(chrom, segment_start_pos, segment_end_pos, read_callback='all')
        df_baseCt = pd.DataFrame({'ref': [n for n in ref_seq], 
                                  'ref_pos': list(range(segment_start_pos, segment_end_pos)), 
                                  'A': A_ct, 'C': C_ct, 'G': G_ct, 'T': T_ct})

        if len(df_baseCt) == 0:
                continue

        df_baseCt['pos_ct'] = 1

        if args.explore_mode:
            # add accessibility score intervals
            df_acc2 = df_acc.loc[(df_acc['end'] > segment_start_pos) & (df_acc['start'] < segment_end_pos)]
            df_acc2['ref_pos'] = df_acc2.apply(lambda row: [x for x in range(row['start'], row['end'])], axis=1)
            df_baseCt = df_baseCt.merge(df_acc2[['ref_pos', 'bin']].explode('ref_pos'), how='left')

            if args.no_ref_n:
                df_baseCt = df_baseCt.melt(id_vars=['pos_ct', 'bin', 'ref'], value_vars=['A', 'C', 'G', 'T'], var_name='alt', value_name='read_ct')
                df_baseCt['base_change'] = df_baseCt['ref'] + '2' + df_baseCt['alt']
                df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['bin', 'base_change', 'pos_ct', 'read_ct']]], axis=0, ignore_index=True)
                df_baseCt_grouped = df_baseCt_grouped.groupby(['bin', 'base_change']).agg({'pos_ct': 'sum', 'read_ct': 'sum'}).reset_index()
                logging.info(f'Writing binned edit data to {args.baseCt_grouped}')
                df_baseCt_grouped.to_csv(args.baseCt_grouped)
            else:
                # get motif for each ref base, DO NOT correct orientation
                df_baseCt['ref_n3'] = [fasta.fetch(chrom, x-motif_span, x+motif_span+1) for x in df_baseCt['ref_pos']]
                df_baseCt = df_baseCt.loc[~df_baseCt['ref_n3'].str.contains('N', regex=False)]

                df_baseCt = df_baseCt.melt(id_vars=['pos_ct', 'bin', 'ref_n3', 'ref'], value_vars=['A', 'C', 'G', 'T'], var_name='alt', value_name='read_ct')
                df_baseCt['base_change'] = df_baseCt['ref'] + '2' + df_baseCt['alt']
                df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['bin', 'ref_n3', 'base_change', 'pos_ct', 'read_ct']]], axis=0, ignore_index=True)
                df_baseCt_grouped = df_baseCt_grouped.groupby(['bin', 'ref_n3', 'base_change']).agg({'pos_ct': 'sum', 'read_ct': 'sum'}).reset_index()
                logging.info(f'Writing binned edit data to {args.baseCt_grouped}')
                df_baseCt_grouped.to_csv(args.baseCt_grouped)

        else:
            df_baseCt['total_ct'] = df_baseCt['A'] + df_baseCt['C'] + df_baseCt['G'] + df_baseCt['T']
            df_baseCt = df_baseCt.loc[(df_baseCt['ref'].isin(['C', 'G'])) & (df_baseCt['total_ct']>0)]                
            df_baseCt['edit_ct'] = df_baseCt.apply(lambda row: row['T'] if row['ref'] == 'C' else row['A'], axis=1)

            # add accessibility score intervals
            df_acc2 = df_acc.loc[(df_acc['end'] > segment_start_pos) & (df_acc['start'] < segment_end_pos)]
            df_acc2['ref_pos'] = df_acc2.apply(lambda row: [x for x in range(row['start'], row['end'])], axis=1)
            df_baseCt = df_baseCt.merge(df_acc2[['ref_pos', 'bin']].explode('ref_pos'), how='left')

            if not args.no_ref_n:
                # get motif for each ref base, correct orientation
                df_baseCt['ref_n3'] = [fasta.fetch(chrom, x-motif_span, x+motif_span+1) for x in df_baseCt['ref_pos']]
                df_baseCt['ref_n3'] = df_baseCt.apply(lambda row: row['ref_n3'] if row['ref'] == 'C' else Seq(row['ref_n3']).reverse_complement().__str__(), axis=1)
                df_baseCt = df_baseCt.loc[~df_baseCt['ref_n3'].str.contains('N', regex=False)]
                
                df_baseCt = df_baseCt[['bin', 'ref_n3', 'pos_ct', 'edit_ct', 'total_ct']].copy()
                df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt], axis=0, ignore_index=True)
                df_baseCt_grouped = df_baseCt_grouped.groupby(['bin', 'ref_n3']).agg({'pos_ct': 'sum', 'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
                logging.info(f'Writing binned edit data to {args.baseCt_grouped}')
                df_baseCt_grouped.to_csv(args.baseCt_grouped)
            else:
                df_baseCt = df_baseCt[['bin', 'pos_ct', 'edit_ct', 'total_ct']].copy()
                df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt], axis=0, ignore_index=True)
                df_baseCt_grouped = df_baseCt_grouped.groupby(['bin']).agg({'pos_ct': 'sum', 'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
                logging.info(f'Writing binned edit data to {args.baseCt_grouped}')
                df_baseCt_grouped.to_csv(args.baseCt_grouped)
 

if not args.explore_mode:
    df_baseCt_grouped['edit_frac'] = df_baseCt_grouped['edit_ct'] / df_baseCt_grouped['total_ct']


if args.baseCt_grouped:
    logging.info(f'Writing binned edit data to {args.baseCt_grouped}')
    df_baseCt_grouped.to_csv(args.baseCt_grouped)























