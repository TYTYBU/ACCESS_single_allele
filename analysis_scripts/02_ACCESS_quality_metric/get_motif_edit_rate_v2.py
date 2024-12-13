import os
import argparse
import pandas as pd
import pysam
from Bio.Seq import Seq
from plotnine import *
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def parse_args():
    parser = argparse.ArgumentParser(
        description='This script calculates Ddd enzyme edit fractions for TF motif surrounding area.',
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

    parser.add_argument('peakDir', 
                        type=str, 
                        default=None,
                        help=('Directory of BED file containing motif peaks. \n'
                              'Default: None'))

    parser.add_argument('--peakFileName_components', 
                        type=str, 
                        default='experiment_id,TF_id,cell_line',
                        help=("Components of BED filename seprated by '_'. \n"
                              "Default: experiment_id,TF_id,cell_line"))

    parser.add_argument('--peakScore_ascending', 
                        type=bool, 
                        default=True,
                        help=('Direction to sort peak scores. \n'
                              "Default: True"))
        
    parser.add_argument('--baseCt_grouped',
                        type=str,
                        default=None,
                        help=('Output composite edit fraction for each motif and cell line as CSV. \n'
                              'Default: None'))
    
    parser.add_argument('--flank_len',
                        type=int,
                        default=500,
                        help=('Length of the flanking region in one direction of the motif. \n'
                              'Default: 500'))

    parser.add_argument('--use_strand',
                        action="store_true", 
                        default=False,
                        help=('Relative position sensitive to motif strand. \n'
                              'Default: False'))

    parser.add_argument('--no_ref_n',
                        action="store_true", 
                        default=False,
                        help=('Do not separate trinucleotide edit motifs. \n'
                              'Default: False'))

    parser.add_argument('--explore_mode',
                        action="store_true", 
                        default=False,
                        help=('Counting all possible base changes. Do not calculate edit fractions. Do not convert orientations for trinucleotide edit motifs. \n'
                              'Default: False'))


    # parser.add_argument('--composite_plot',
    #                     type=str,
    #                     default=None,
    #                     help=('Draw composite edit fraction plot for motifs. \n'
    #                           'Default: None')) 

    # parser.add_argument('--plot_height',
    #                     type=int,
    #                     default=6,
    #                     help=('Height of composite plot. \n'
    #                           'Default: 6'))

    # parser.add_argument('--plot_width',
    #                     type=int,
    #                     default=12,
    #                     help=('Width of composite plot. \n'
    #                           'Default: 12'))

    # parser.add_argument('--bamSuffix',
    #                     type=str,
    #                     default='',
    #                     help=("BAM file suffix to remove from sample name. \n"
    #                           "Default: ''"))

    return parser.parse_args()


# load arguments
args = parse_args()
bam_path=args.bamFile
fasta_path=args.fastaFile
dir_path = args.peakDir
# sample_id = os.path.basename(bam_path).replace(args.bamSuffix, '')
flank_len = args.flank_len

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
if not os.path.exists(fasta_path + '.fai'):
    print("fasta file need to indexed first!")
    exit(1)
fasta = pysam.FastaFile(fasta_path)

# load factorbook TF motif peaks
logging.info(f'Loading TF motif peaks from {dir_path}')
dir_list = os.listdir(dir_path)
dir_list = [x for x in dir_list if ('.txt' in x) or ('.bed' in x)]
df = pd.DataFrame([x.replace('.gz', '').replace('.txt', '').replace('.bed', '').split('_') for x in dir_list], columns=args.peakFileName_components.split(','))
df['file_name'] = dir_list
if 'cell_line' not in df.columns:
    df['cell_line'] = 'unknown'

# compile all TF motif locations
logging.info(f'Compiling TF motif regions with flanking length: {flank_len}nt.')
if args.peakScore_ascending:
    logging.info(f'TF motifs are sorted by ascending scores.')
else:
    logging.info(f'TF motifs are sorted by descending scores.')

major_chrs = ['chr'+y for y in ([str(x) for x in range(1, 23)] + ['X', 'Y'])]

df_TF_all = pd.DataFrame()
for index, row in df.iterrows():
    in_file_path = '/'.join([dir_path, row['file_name']])
    df_TF = pd.read_csv(in_file_path, sep='\t', header=None)
    df_TF = df_TF.iloc[:, :5]
    df_TF.columns = ["chr", "start", "end", "strand", "score"]

    df_TF['TF_name'] = df_TF.apply(lambda row: '_'.join([row['chr'], str(row['start']), str(row['end'])]), axis=1)
    df_TF['TF_center'] = df_TF.apply(lambda row: int((row['start'] + row['end'])/2), axis=1)
    df_TF['flank_start'] = df_TF['TF_center'] - flank_len
    df_TF['flank_end'] = df_TF['TF_center'] + flank_len
    df_TF['TF_id'] = row['TF_id']
    df_TF['cell_line'] = row['cell_line']

    # currently the pipeline only allow TF motifs on major human chromosomes
    df_TF = df_TF.loc[df_TF['chr'].isin(major_chrs)]

    df_TF = df_TF.sort_values(by=['score'], ignore_index=True, ascending=args.peakScore_ascending).head(10000)
    df_TF_all = pd.concat([df_TF_all, df_TF], axis=0, ignore_index=True)

logging.info(f'Total motif region created: {df_TF_all.shape[0]}')


# calculate edit frequency for each 3mer across motifs
logging.info(f'Calculating edit frequency across all motifs regions ...')
ct=0
df_baseCt_grouped = pd.DataFrame()
for index, row in df_TF_all.iterrows():
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']
    TF_center = row['TF_center']
    strand = row['strand']
    
    ref_seq = fasta.fetch(chr, flank_start, flank_end)
    (A_ct, C_ct, G_ct, T_ct) = bam.count_coverage(chr, flank_start, flank_end, read_callback='nofilter')
    
    df_baseCt = pd.DataFrame({'ref': [n for n in ref_seq], 
                              'ref_pos': list(range(flank_start,flank_end)), 
                              'A': A_ct, 'C': C_ct, 'G': G_ct, 'T': T_ct})

    if len(df_baseCt) == 0:
            continue

    if args.use_strand:
        if strand == "+":
            df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center
        else:
            df_baseCt['relative_pos'] = TF_center - df_baseCt['ref_pos']
    else:
        df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center

    if args.explore_mode:
        if args.no_ref_n:
            df_baseCt = df_baseCt.melt(id_vars=['relative_pos', 'ref'], value_vars=['A', 'C', 'G', 'T'], var_name='alt', value_name='read_ct')
            df_baseCt['base_change'] = df_baseCt['ref'] + '2' + df_baseCt['alt']
            df_baseCt['TF_id'] = row['TF_id']
            df_baseCt['cell_line'] = row['cell_line']
            df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['TF_id', 'cell_line', 'relative_pos', 'base_change', 'read_ct']]], axis=0, ignore_index=True)
        else:
            # get motif for each ref base, correct orientation
            df_baseCt['ref_n3'] = [fasta.fetch(chr, x-1, x+2) for x in df_baseCt['ref_pos']]
            df_baseCt = df_baseCt.loc[~df_baseCt['ref_n3'].str.contains('N', regex=False)]

            df_baseCt = df_baseCt.melt(id_vars=['relative_pos', 'ref_n3', 'ref'], value_vars=['A', 'C', 'G', 'T'], var_name='alt', value_name='read_ct')
            df_baseCt['base_change'] = df_baseCt['ref'] + '2' + df_baseCt['alt']
            df_baseCt['TF_id'] = row['TF_id']
            df_baseCt['cell_line'] = row['cell_line']
            df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['TF_id', 'cell_line', 'relative_pos', 'ref_n3', 'base_change', 'read_ct']]], axis=0, ignore_index=True)

        ct += 1
        if ct % 1000 == 0:
            if not args.no_ref_n:
                df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'base_change', 'ref_n3']).agg({'read_ct': 'sum'}).reset_index()
            else:
                df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'base_change']).agg({'read_ct': 'sum'}).reset_index()
            logging.info(f'Processed {ct} motifs ...')
    else:
        df_baseCt['total_ct'] = df_baseCt['A'] + df_baseCt['C'] + df_baseCt['G'] + df_baseCt['T']
        df_baseCt = df_baseCt.loc[df_baseCt['ref'].isin(['C', 'G'])]
        df_baseCt['edit_ct'] = df_baseCt.apply(lambda row: row['T'] if row['ref'] == 'C' else row['A'], axis=1)
        
        if not args.no_ref_n:
            # get 3mer for each ref base, correct orientation
            df_baseCt['ref_n3'] = [fasta.fetch(chr, x-1, x+2) for x in df_baseCt['ref_pos']]
            df_baseCt['ref_n3'] = df_baseCt.apply(lambda row: row['ref_n3'] if row['ref'] == 'C' else Seq(row['ref_n3']).reverse_complement().__str__(), axis=1)
            df_baseCt = df_baseCt.loc[~df_baseCt['ref_n3'].str.contains('N', regex=False), ['relative_pos', 'ref_n3', 'edit_ct', 'total_ct']]

        df_baseCt['TF_id'] = row['TF_id']
        df_baseCt['cell_line'] = row['cell_line']
        df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt], axis=0, ignore_index=True)

        ct += 1
        if ct % 1000 == 0:
            if not args.no_ref_n:
                df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'ref_n3']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
            else:
                df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
            logging.info(f'Processed {ct} motifs ...')

if args.explore_mode:
    if not args.no_ref_n:
        df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'base_change', 'ref_n3']).agg({'read_ct': 'sum'}).reset_index()
    else:
        df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'base_change']).agg({'read_ct': 'sum'}).reset_index()

else:
    if not args.no_ref_n:
        df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos', 'ref_n3']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
    else:
        df_baseCt_grouped = df_baseCt_grouped.groupby(['TF_id', 'cell_line', 'relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()

    df_baseCt_grouped = df_baseCt_grouped.loc[df_baseCt_grouped['total_ct']>0].reset_index(drop=True)
    df_baseCt_grouped['edit_frac'] = df_baseCt_grouped['edit_ct'] / df_baseCt_grouped['total_ct']


if args.baseCt_grouped:
    logging.info(f'Writing motif edit data to {args.baseCt_grouped}')
    df_baseCt_grouped.to_csv(args.baseCt_grouped, index=False)



# # plot edit fraction for each TF and cell line
# if args.composite_plot:
#     logging.info(f'Generating composite edit rate plot at {args.composite_plot}')
#     plt_title = sample_id
#     plt_subtitle = "Ddd1 TFs 10K-rankset composite edit fraction"

#     if not args.no_ref_n:
#         plt = (
#             ggplot(df_baseCt_grouped, aes(x='relative_pos', y='edit_frac', color='ref_n3'))
#             + geom_point(size=0.3, alpha=0.2)
#             + stat_smooth(method='mavg', method_args={'window': 5}, size=0.25, se=False)
#             + scale_color_manual(values=['grey']*12 + ['green', 'blue', 'orange', 'red'], name='3mer')
#             + labs(x='Relative position', y='Edit fraction', title=plt_title, subtitle=plt_subtitle)
#             + facet_grid('TF_id ~ cell_line', scales='free_y')
#             + theme_light()
#             + theme(legend_key=element_rect(color = "white"))
#         )
#         plt.save(filename=args.composite_plot, width=args.plot_width, height=args.plot_height)
#     else:
#         plt = (
#             ggplot(df_baseCt_grouped, aes(x='relative_pos', y='edit_frac'))
#             + geom_point(size=0.3, alpha=0.1)
#             + stat_smooth(method='mavg', method_args={'window': 6}, size=0.25, se=False)
#             + labs(x='Relative position', y='Edit fraction', title=plt_title, subtitle=plt_subtitle)
#             + facet_grid('TF_id ~ cell_line', scales='free_y')
#             + theme_light()
#             + theme(legend_key=element_rect(color = "white"))
#         )
#         plt.save(filename=args.composite_plot, width=args.plot_width, height=args.plot_height)
    













