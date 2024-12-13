import pandas as pd
import numpy as np
from Bio.Seq import Seq
import pysam
from tqdm import tqdm
import tensorflow as tf
from scipy.stats import spearmanr, pearsonr
import ast
from scipy.fftpack import fft

def get_edit_pos(aligned_pairs = None,
                 q_seq: str = None,
                 all_pos: bool = True):
    edit_pos = tuple()
    edit_type = tuple()
    r_start = None
    r_end = None
    
    for q_pos, r_pos, r_base in aligned_pairs:
        r_pos = int(r_pos)
        r_base = r_base.upper()
        q_base = q_seq[q_pos]

        if r_base == 'C' and q_base == 'T':
            edit_pos = edit_pos + (r_pos,)
            edit_type = edit_type + ('C2T',)
        elif all_pos and r_base == 'C' and q_base == 'C':
            edit_pos = edit_pos + (r_pos,)
            edit_type = edit_type + ('C2C',)
        elif r_base == 'G' and q_base == 'A':
            edit_pos = edit_pos + (r_pos,)
            edit_type = edit_type + ('G2A',)
        elif all_pos and r_base == 'G' and q_base == 'G':
            edit_pos = edit_pos + (r_pos,)
            edit_type = edit_type + ('G2G',)
    return edit_pos, edit_type

def calc_read_level(df_read: pd.DataFrame = None):
    temp = [0]*len(df_read)
    read_level = [0]*len(df_read)
    for i in range(len(df_read)):
        if i>0:
            while df_read['R1__r_start'][i] < temp[read_level[i]]:
                read_level[i]+=1
        temp[read_level[i]] = max(df_read['R1__r_end'][i], df_read['R2__r_end'][i])

    return read_level


def get_read_df(bam = None,
                chr: str = None, 
                start: int = None, 
                end: int = None,
                TF_start: int = None,
                read_level: bool = True) -> pd.DataFrame:
    df_read = pd.DataFrame(columns = ['read_id', 'r_start', 'r_end', 'edit_pos', 'edit_type'])
    i=0
    for read in bam.fetch(chr, start, end):
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
        q_seq = read.query_sequence.upper()
        temp = read.query_name.split(':')
        read_id = ''.join(temp[4:7])
        all_r_pos = read.get_reference_positions()
        r_start = min(all_r_pos)
        r_end = max(all_r_pos)

        # get relative edit positions
        edit_pos, edit_type = get_edit_pos(aligned_pairs, q_seq)
        edit_pos = tuple(x-TF_start for x in edit_pos)
        
        df_read.loc[i] = [read_id, r_start, r_end, edit_pos, edit_type]
        i+=1

    # separate paired reads and singleton reads
    col_names = df_read.columns
    col_names_singleton = [col_names[0]] + ['R1__' + x for x in col_names[1:]]
    col_names_paired = [col_names[0]] + ['R1__' + x for x in col_names[1:]] + ['R2__' + x for x in col_names[1:]] 
    
    df_singleton = pd.DataFrame(columns = col_names_singleton)
    df_paired = pd.DataFrame(columns = col_names_paired)
    
    i=0
    for read_id in set(df_read['read_id']):
        temp = df_read.loc[df_read['read_id'] == read_id].reset_index(drop=True)
        if (len(temp) == 1):
            df_singleton.loc[i] = temp.loc[0].tolist()
        else:
            temp1 = temp.loc[0].tolist()
            temp2 = temp.loc[1].tolist()
            if temp2[1] < temp1[1]:
                temp3 = temp1
                temp1 = temp2
                temp2 = temp3
            df_paired.loc[i] = temp1 + temp2[1:]
        i += 1

    # merge paired reads and singleton reads
    df_read = df_paired.merge(df_singleton, how='outer')
    df_read.sort_values(by=['R1__r_start'], inplace=True, ignore_index=True)
    
    # calculate the position of read on y axis
    if read_level:
        df_read['read_level'] = calc_read_level(df_read)
    
    # relative position to TF start
    df_read['R1__r_start'] = df_read['R1__r_start'] - TF_start
    df_read['R1__r_end'] = df_read['R1__r_end'] - TF_start
    df_read['R2__r_start'] = df_read['R2__r_start'] - TF_start
    df_read['R2__r_end'] = df_read['R2__r_end'] - TF_start

    return df_read

def get_edits_df(df_read: pd.DataFrame = None) -> pd.DataFrame:
    df_edits = pd.DataFrame(columns = ['x_pos', 'y_pos', 'group'])
    j=0
    for index, row in df_read.iterrows():
        R1__edit_pos = row['R1__edit_pos']
        R1__edit_type = row['R1__edit_type']
        R2__edit_pos = row['R2__edit_pos']
        R2__edit_type = row['R2__edit_type']
        read_level = row['read_level']

        if isinstance(R1__edit_pos, tuple) and len(R1__edit_pos)>0:
            df_temp = pd.DataFrame({'x_pos': list(R1__edit_pos), 'y_pos': read_level, 'group': R1__edit_type})
            df_edits = pd.concat([df_edits, df_temp], ignore_index=True)

        if isinstance(R2__edit_pos, tuple) and len(R2__edit_pos)>0:
            df_temp = pd.DataFrame({'x_pos': list(R2__edit_pos), 'y_pos': read_level, 'group': R2__edit_type})
            df_edits = pd.concat([df_edits, df_temp], ignore_index=True)

    df_edits['x_pos'] = df_edits['x_pos'].to_list()
    df_edits['y_pos'] = df_edits['y_pos'].to_list()
    df_edits['group'] = pd.Categorical(df_edits['group'], categories=['C2T', 'G2A', 'C2C', 'G2G'], ordered=True)
                
    return df_edits




def get_flank_pos(df, flank_len=300):
    df['flank_start'] = df['start'] - flank_len
    df['flank_end'] = df['end'] + flank_len
    df['center'] = df.apply(lambda row: int((row['start']+row['end'])/2) if row['strand']=='+' else int((row['start']+row['end'])/2+0.5), axis=1)
    return(df)

def df_motifs_2_composite_edit_frac(bam, fasta, df_motifs, editable_pos_only=False, use_strand=True):
    if not all([x in df_motifs.columns for x in ['flank_start', 'flank_end', 'center']]):
        df_motifs = get_flank_pos(df_motifs)

    ct=0
    df_baseCt_grouped = pd.DataFrame()
    for index, row in df_motifs.iterrows():
        chr = row['chr']
        flank_start = row['flank_start']
        flank_end = row['flank_end']
        TF_center = row['center']
        strand = row['strand']
        
        ref_seq = fasta.fetch(chr, flank_start, flank_end)
        (A_ct, C_ct, G_ct, T_ct) = bam.count_coverage(chr, flank_start, flank_end, read_callback='nofilter')
        
        df_baseCt = pd.DataFrame({'ref_base': [n for n in ref_seq], 
                                'ref_pos': list(range(flank_start,flank_end)), 
                                'A_ct': A_ct, 'C_ct': C_ct, 'G_ct': G_ct, 'T_ct': T_ct})
        if use_strand:
            if strand == "+":
                df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center
            else:
                df_baseCt['relative_pos'] = TF_center - df_baseCt['ref_pos']
        else:
            df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center
            
        df_baseCt = df_baseCt.loc[df_baseCt['ref_base'].isin(['C', 'G'])]
        df_baseCt['edit_ct'] = df_baseCt.apply(lambda row: row['T_ct'] if row['ref_base'] == 'C' else row['A_ct'], axis=1)
        if editable_pos_only:
            df_baseCt['total_ct'] = df_baseCt.apply(lambda row: (row['T_ct']+row['C_ct']) if row['ref_base'] == 'C' else (row['A_ct']+row['G_ct']), axis=1)
        else:
            df_baseCt['total_ct'] = df_baseCt['A_ct'] + df_baseCt['C_ct'] + df_baseCt['G_ct'] + df_baseCt['T_ct']

        df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['relative_pos', 'edit_ct', 'total_ct']]], axis=0, ignore_index=True)

        ct += 1
        if ct % 1000 == 0:
            df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()

    df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'})
    df_baseCt_grouped = df_baseCt_grouped.reset_index()
    df_baseCt_grouped['edit_frac'] = df_baseCt_grouped['edit_ct'] / df_baseCt_grouped['total_ct']
    return(df_baseCt_grouped)

def motif_row_2_edit_frac(bam, fasta, row):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']
    TF_center = row['center']
    strand = row['strand']
    
    ref_seq = fasta.fetch(chr, flank_start, flank_end)
    (A_ct, C_ct, G_ct, T_ct) = bam.count_coverage(chr, flank_start, flank_end, read_callback='nofilter')
    
    df_baseCt = pd.DataFrame({'ref_base': [n for n in ref_seq], 
                            'ref_pos': list(range(flank_start,flank_end)), 
                            'A_ct': A_ct, 'C_ct': C_ct, 'G_ct': G_ct, 'T_ct': T_ct})
    if strand == "+":
        df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center
    else:
        df_baseCt['relative_pos'] = TF_center - df_baseCt['ref_pos']
    df_baseCt = df_baseCt.loc[df_baseCt['ref_base'].isin(['C', 'G'])]
    df_baseCt['edit_ct'] = df_baseCt.apply(lambda row: row['T_ct'] if row['ref_base'] == 'C' else row['A_ct'], axis=1)
    df_baseCt['total_ct'] = df_baseCt['A_ct'] + df_baseCt['C_ct'] + df_baseCt['G_ct'] + df_baseCt['T_ct']
    df_baseCt['edit_frac'] = df_baseCt['edit_ct'] / df_baseCt['total_ct']
    return(df_baseCt[['relative_pos', 'edit_ct', 'total_ct', 'edit_frac']])

def df_motifs_2_edit_PFM(bam, fasta, df_motifs):
    if not all([x in df_motifs.columns for x in ['flank_start', 'flank_end', 'center']]):
        df_motifs = get_flank_pos(df_motifs)

    ct=0
    df_baseCt_grouped = pd.DataFrame()
    for index, row in df_motifs.iterrows():
        chr = row['chr']
        flank_start = row['flank_start']
        flank_end = row['flank_end']
        TF_center = row['center']
        strand = row['strand']
        
        ref_seq = fasta.fetch(chr, flank_start, flank_end)
        (A_ct, C_ct, G_ct, T_ct) = bam.count_coverage(chr, flank_start, flank_end, read_callback='nofilter')
        
        df_baseCt = pd.DataFrame({'ref_base': [n for n in ref_seq], 
                                'ref_pos': list(range(flank_start,flank_end)), 
                                'A': A_ct, 'C': C_ct, 'G': G_ct, 'T': T_ct})
        if strand == "+":
            df_baseCt['relative_pos'] = df_baseCt['ref_pos'] - TF_center
        else:
            df_baseCt['relative_pos'] = TF_center - df_baseCt['ref_pos']

        # Y = C or T
        df_baseCt['Y'] = df_baseCt.apply(lambda row: row['T'] if row['ref_base'] == 'C' else 0, axis=1)
        df_baseCt['T'] = df_baseCt.apply(lambda row: 0 if row['ref_base'] == 'C' else row['T'], axis=1)

        # R = G or A
        df_baseCt['R'] = df_baseCt.apply(lambda row: row['A'] if row['ref_base'] == 'G' else 0, axis=1)
        df_baseCt['A'] = df_baseCt.apply(lambda row: 0 if row['ref_base'] == 'G' else row['A'], axis=1)

        df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['relative_pos', 'A', 'C', 'G', 'T', 'R', 'Y']]], axis=0, ignore_index=True)

        ct += 1
        if ct % 1000 == 0:
            df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'A': 'sum', 'C': 'sum', 'G': 'sum', 'T': 'sum', 'R': 'sum', 'Y': 'sum'})
            df_baseCt_grouped = df_baseCt_grouped.reset_index()

    df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'A': 'sum', 'C': 'sum', 'G': 'sum', 'T': 'sum', 'R': 'sum', 'Y': 'sum'})
    df_baseCt_grouped = df_baseCt_grouped.reset_index()
    return(df_baseCt_grouped)

def read2edits(read, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A'], penalize_non_edit=False):
    edit_types = [x.split('2') for x in edit_types]
    flank_start = TF_center - flank_len
    flank_end = TF_center + flank_len

    # for a given read, mark edited bases as 1 and non-edited bases as 0
    # output as a list for the flanking region
    aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
    q_seq = read.query_sequence.upper()
    edits = [0] * (flank_len*2+1)

    for q_pos, r_pos, ref in aligned_pairs:
        ref = ref.upper()
        if motif_strand == '+':
            i = r_pos - flank_start
        else:
            i = flank_end - r_pos

        if (i >= 0) and (i <= flank_len*2):
            alt = q_seq[q_pos]
            
            for ref2,alt2 in edit_types:
                if ref == ref2 and alt == alt2:
                    edits[i] = 1
                    break
                if penalize_non_edit and ref == alt == ref2:
                    edits[i] = -1
                    break
    return(edits)

def read2covs(read, TF_center, flank_len, motif_strand, editable_pos_only=False, editable_bases = ['C', 'G']):
    flank_start = TF_center - flank_len
    flank_end = TF_center + flank_len

    # for the flanking region, mark read positions as 1 and non-read positions 0
    covs = [0] * (flank_len*2+1)
    if editable_pos_only:
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
        for q_pos, r_pos, ref in aligned_pairs:
            ref = ref.upper()
            if motif_strand == '+':
                i = r_pos - flank_start
            else:
                i = flank_end - r_pos
            if (i >= 0) and (i <= flank_len*2) and (ref in editable_bases):
                covs[i] = 1
    else:
        if motif_strand == '+':
            i_start = read.reference_start - flank_start
            i_end = read.reference_end - flank_start
        else:
            i_start = flank_end - read.reference_end
            i_end = flank_end - read.reference_start
        for i in range(i_start, i_end):
            if (i >= 0) and (i <= flank_len*2):
                covs[i] = 1
    return(covs)

def motif_row_2_mat_edits(bam, row, flank_len, cover_center=False, penalize_non_edit=False):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']
    TF_center = row['center']
    strand = row['strand']

    mat_edits = []
    for read in bam.fetch(chr, flank_start, flank_end):
        if cover_center:
            if read.reference_start <= TF_center < read.reference_end:
                edits = read2edits(read, TF_center, flank_len, strand, penalize_non_edit)
                mat_edits.append(edits)                
        else:
            edits = read2edits(read, TF_center, flank_len, strand, penalize_non_edit)  
            mat_edits.append(edits)           
    return(np.asarray(mat_edits))

def motif_row_2_mat_covs(bam, row, flank_len, cover_center=False):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']
    TF_center = row['center']
    strand = row['strand']

    mat_covs = []
    for read in bam.fetch(chr, flank_start, flank_end):
        if cover_center:
            if read.reference_start <= TF_center < read.reference_end:
                covs = read2covs(read, TF_center, flank_len, strand)
                mat_covs.append(covs)                
        else:
            covs = read2covs(read, TF_center, flank_len, strand)  
            mat_covs.append(covs)           
    return(np.asarray(mat_covs))

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

def df_motifs_2_mat_edits(bam, df_motifs, flank_len, cover_center=False, penalize_non_edit=False):
    if not all([x in df_motifs.columns for x in ['flank_start', 'flank_end']]):
        df_motifs = get_flank_pos(df_motifs, flank_len)

    mat_edits = []
    for index, row in df_motifs.iterrows():
        chr = row['chr']
        flank_start = row['flank_start']
        flank_end = row['flank_end']
        TF_center = row['center']
        strand = row['strand']

        for read in bam.fetch(chr, flank_start, flank_end):
            if cover_center:
                if read.reference_start <= TF_center < read.reference_end:
                    edits = read2edits(read, TF_center, flank_len, strand, penalize_non_edit)
                    mat_edits.append(edits)
            else:
                edits = read2edits(read, TF_center, flank_len, strand, penalize_non_edit)
                mat_edits.append(edits)        
    return(np.asarray(mat_edits))

def df_motifs_2_mat_covs(bam, df_motifs, flank_len, cover_center=False):
    if not all([x in df_motifs.columns for x in ['flank_start', 'flank_end']]):
        df_motifs = get_flank_pos(df_motifs, flank_len)

    mat_covs = []
    for index, row in df_motifs.iterrows():
        chr = row['chr']
        flank_start = row['flank_start']
        flank_end = row['flank_end']
        TF_center = row['center']
        strand = row['strand']

        for read in bam.fetch(chr, flank_start, flank_end):
            if cover_center:
                if read.reference_start <= TF_center < read.reference_end:
                    covs = read2covs(read, TF_center, flank_len, strand)
                    mat_covs.append(covs)
            else:
                covs = read2covs(read, TF_center, flank_len, strand)
                mat_covs.append(covs)        
    return(np.asarray(mat_covs))

def merge_lists(list_a, list_b):
    merged_list = [0] * len(list_a)
    for i in range(len(list_a)):
        if list_a[i] == list_b[i]:
            merged_list[i] = list_a[i]
        elif (list_a[i] == 1 and list_b[i] == 0) or (list_a[i] == 0 and list_b[i] == 1):
            merged_list[i] = 1
        elif (list_a[i] == -1 and list_b[i] == 0) or (list_a[i] == 0 and list_b[i] == -1):
            merged_list[i] = -1
        else:
            merged_list[i] = 0
    return merged_list
    
def merge_read_pair_edits(read1, read2, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A'], penalize_non_edit=False):
    if type(read1) != pysam.AlignedSegment:
        return read2edits(read2, TF_center, flank_len, motif_strand, edit_types, penalize_non_edit)
    elif type(read2) != pysam.AlignedSegment:
        return read2edits(read1, TF_center, flank_len, motif_strand, edit_types, penalize_non_edit)
    else:
        edits1 = read2edits(read1, TF_center, flank_len, motif_strand, edit_types, penalize_non_edit)
        edits2 = read2edits(read2, TF_center, flank_len, motif_strand, edit_types, penalize_non_edit)
        return(merge_lists(edits1, edits2))
    
def merge_read_pair_covs(read1, read2, TF_center, flank_len, motif_strand, editable_pos_only=False, editable_bases = ['C', 'G']):
    if type(read1) != pysam.AlignedSegment:
        return read2covs(read2, TF_center, flank_len, motif_strand, editable_pos_only, editable_bases)
    elif type(read2) != pysam.AlignedSegment:
        return read2covs(read1, TF_center, flank_len, motif_strand, editable_pos_only, editable_bases)
    else:
        covs1 = read2covs(read1, TF_center, flank_len, motif_strand, editable_pos_only, editable_bases)
        covs2 = read2covs(read2, TF_center, flank_len, motif_strand, editable_pos_only, editable_bases)
        return(merge_lists(covs1, covs2))

def motif_row_2_df_read_pairs(bam, row):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']

    r1_names = []
    r1_reads = []
    r2_names = []
    r2_reads = []
    for read in bam.fetch(chr, flank_start, flank_end):
        read_name = read.query_name
        if read.is_read2:
            r2_names.append(read_name)
            r2_reads.append(read)
        else:
            r1_names.append(read_name)
            r1_reads.append(read)
    
    df_r1 = pd.DataFrame({'name': r1_names, 'r1_reads': r1_reads})
    df_r2 = pd.DataFrame({'name': r2_names, 'r2_reads': r2_reads})
    df_read_pairs = pd.merge(df_r1, df_r2, on='name', how='outer')
    return(df_read_pairs)

def edits_frac_in_range(edits, covs, relative_pos_list, flank_len, verbose=False):
    idx_list = [i+flank_len for i in relative_pos_list]
    edits_in_range = [edits[i] for i in idx_list if covs[i]]
    if verbose:
        print(edits_in_range)
    if len(edits_in_range):
        return(np.mean(edits_in_range))
    else:
        return(float('nan'))

def covs_frac_in_range(covs, relative_pos_list, flank_len, verbose=False):
    idx_list = [i+flank_len for i in relative_pos_list]
    covs_in_range = [covs[i] for i in idx_list]
    if verbose:
        print(covs_in_range)
    if len(covs_in_range):
        return(np.mean(covs_in_range))
    else:
        return(float('nan'))

def df_read_pairs_2_composite_edit_frac(df_read_pairs, flank_len, edit_types = ['C2T', 'G2A']):
    edit_types = [x.split('2') for x in edit_types]
    edit_ct = [0] * (flank_len*2+1)
    total_ct = [0] * (flank_len*2+1)

    for index, row in df_read_pairs.iterrows():
        TF_center = row['center']
        strand = row['strand']
        flank_start = TF_center - flank_len
        flank_end = TF_center + flank_len

        for read in row[['r1_reads', 'r2_reads']]:
            if type(read) == pysam.AlignedSegment:
                aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
                q_seq = read.query_sequence.upper()

                for q_pos, r_pos, ref in aligned_pairs:
                    ref = ref.upper()
                    if strand == '+':
                        i = r_pos - flank_start
                    else:
                        i = flank_end - r_pos

                    if (i >= 0) and (i <= flank_len*2):
                        alt = q_seq[q_pos]
                        for ref2,alt2 in edit_types:
                            if ref == ref2:
                                total_ct[i] += 1
                                if alt == alt2:
                                    edit_ct[i] += 1
                                break

    df_baseCt = pd.DataFrame({'relative_pos': [i for i in range(-flank_len, flank_len+1)],
                              'edit_ct': edit_ct,
                              'total_ct': total_ct})
    df_baseCt['edit_frac'] = df_baseCt['edit_ct'] / df_baseCt['total_ct']
    return(df_baseCt)

def mat_edits_2_edit_frac(mat_edits, penalize_non_edit=False):
    flank_len = int((mat_edits.shape[1] - 1)/2)
    edit_frac = [float('nan')] * mat_edits.shape[1]
    unedit_value = -1 if penalize_non_edit else 0

    for i in range(mat_edits.shape[1]):
        uniq_value, value_ct = np.unique(mat_edits[:,i], return_counts=True)
        edit_ct = value_ct[uniq_value == 1]
        edit_ct = edit_ct[0] if len(edit_ct) else 0
        unedit_ct = value_ct[uniq_value == unedit_value]
        unedit_ct = unedit_ct[0] if len(unedit_ct) else 0
        edit_frac[i] = edit_ct / (edit_ct + unedit_ct)
    
    df_edit_frac = pd.DataFrame({'relative_pos': [x for x in range(-flank_len, flank_len+1)], 'edit_frac': edit_frac})
    return(df_edit_frac)

def mat_edits_2_edit_frac_by_read_type_labels(mat_edits, read_type_labels, penalize_non_edit=False, read_type_list = ['unbound', 'bound', 'recently bound']):
    df_composite_edit_frac = pd.DataFrame()
    for i in range(len(read_type_list)):
        mat_edits_label = mat_edits[read_type_labels == i,:]
        df_temp = mat_edits_2_edit_frac(mat_edits_label, penalize_non_edit)
        df_temp['read_type'] = read_type_list[i]
        df_composite_edit_frac = pd.concat([df_composite_edit_frac, df_temp])
    return(df_composite_edit_frac)

def df_motifs_2_read_pairs_2_mat_edits(bam, df_motifs, flank_len=150, penalize_non_edit=False, 
                                       center_trough_pos_list = [x for x in range(-10,11)],
                                       left_peak_pos_list = [x for x in range(-75,-9)],
                                       right_peak_pos_list = [x for x in range(10,76)],
                                       center_trough_covs_frac_thres = 0.5,
                                       peak_covs_frac_thres = 0.5):
    mat_edits = []
    for index, row in df_motifs.iterrows():
        TF_center = row['center']
        strand = row['strand']
        df_read_pairs = motif_row_2_df_read_pairs(bam, row)

        if len(df_read_pairs):
            df_read_pairs['chr'] = row['chr']
            df_read_pairs['center'] = TF_center
            df_read_pairs['strand'] = strand

            # get overall coverage list, filter reads by coverage in feature regions
            df_read_pairs['covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
            df_read_pairs['center_trough_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
            df_read_pairs['left_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, left_peak_pos_list, flank_len))
            df_read_pairs['right_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, right_peak_pos_list, flank_len))
            df_read_pairs = df_read_pairs.loc[(df_read_pairs['center_trough_covs_frac']>=center_trough_covs_frac_thres) & 
                                            ((df_read_pairs['left_peak_covs_frac']>=peak_covs_frac_thres) | (df_read_pairs['right_peak_covs_frac']>=peak_covs_frac_thres))]

            if len(df_read_pairs):
                df_read_pairs['edits_all'] = df_read_pairs.apply(lambda row2: merge_read_pair_edits(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand, penalize_non_edit=penalize_non_edit), axis=1)
                mat_edits += df_read_pairs['edits_all'].tolist()

    return(np.asarray(mat_edits))

def motif_row_2_read_pairs_2_mat_edits(bam, row, flank_len=150, penalize_non_edit=False, 
                                       center_trough_pos_list = [x for x in range(-10,11)],
                                       left_peak_pos_list = [x for x in range(-75,-9)],
                                       right_peak_pos_list = [x for x in range(10,76)],
                                       center_trough_covs_frac_thres = 0.5,
                                       peak_covs_frac_thres = 0.5):
    mat_edits = []
    TF_center = row['center']
    strand = row['strand']
    df_read_pairs = motif_row_2_df_read_pairs(bam, row)

    if len(df_read_pairs):
        df_read_pairs['chr'] = row['chr']
        df_read_pairs['center'] = TF_center
        df_read_pairs['strand'] = strand

        # get overall coverage list, filter reads by coverage in feature regions
        df_read_pairs['covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
        df_read_pairs['center_trough_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
        df_read_pairs['left_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, left_peak_pos_list, flank_len))
        df_read_pairs['right_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, right_peak_pos_list, flank_len))
        df_read_pairs = df_read_pairs.loc[(df_read_pairs['center_trough_covs_frac']>=center_trough_covs_frac_thres) & 
                                        ((df_read_pairs['left_peak_covs_frac']>=peak_covs_frac_thres) | (df_read_pairs['right_peak_covs_frac']>=peak_covs_frac_thres))]

        if len(df_read_pairs):
            df_read_pairs['edits_all'] = df_read_pairs.apply(lambda row2: merge_read_pair_edits(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand, penalize_non_edit=penalize_non_edit), axis=1)
            mat_edits += df_read_pairs['edits_all'].tolist()

    return(np.asarray(mat_edits))

def mat_edits_1d_to_2d(mat_edits, progress_bar=False):
    # function to convert the 1d edit array [edit, non-edit, non-editable] to 2d edit array [[edits],[editables]]
    # assumes that non-edits at editable positions are penalized in mat_edits_1d
    edits_2d = []
    if progress_bar:
        for edits_1d in mat_edits:
            covs = [(1 if x != 0 else 0) for x in edits_1d]
            edits = [(1 if x == 1 else 0) for x in edits_1d]
            edits_2d.append([edits, covs])
    else:
        for edits_1d in mat_edits:
            covs = [(1 if x != 0 else 0) for x in edits_1d]
            edits = [(1 if x == 1 else 0) for x in edits_1d]
            edits_2d.append([edits, covs])
    return(np.asarray(edits_2d))

def mat_edits_2d_2_edit_frac(mat_edits):
    mat_edits = mat_edits.squeeze(axis=3)
    flank_len = int((mat_edits.shape[2] - 1)/2)
    edit_frac = [float('nan')] * mat_edits.shape[2]

    for i in range(mat_edits.shape[2]):
        edits = mat_edits[:,0,i]
        covs = mat_edits[:,1,i]
        edit_frac[i] = np.sum(edits) / np.sum(covs)
    
    df_edit_frac = pd.DataFrame({'relative_pos': [x for x in range(-flank_len, flank_len+1)], 'edit_frac': edit_frac})
    return(df_edit_frac)

def mat_edits_2d_2_edit_frac_by_read_type_labels(mat_edits, read_type_labels, read_type_list = ['unbound', 'bound', 'recently bound']):
    df_composite_edit_frac = pd.DataFrame()
    for i in range(len(read_type_list)):
        mat_edits_label = mat_edits[read_type_labels == i,:]
        df_temp = mat_edits_2d_2_edit_frac(mat_edits_label)
        df_temp['read_type'] = read_type_list[i]
        df_composite_edit_frac = pd.concat([df_composite_edit_frac, df_temp])
    return(df_composite_edit_frac)

def df_motifs_bin_by_score(df_motifs, score_col='chipseq_score', breaks=[x for x in range(0,401,50)]+[float('inf')], sample_size=1000):
    # assumes df_motifs have motif_type column
    df_inactive = df_motifs.loc[df_motifs['motif_type'] == 'Inactive']
    df_active = df_motifs.loc[df_motifs['motif_type'] == 'Active']

    score_bins = ['Inactive']
    df_motifs_bin = df_inactive.sample(sample_size) if len(df_inactive)>sample_size else df_inactive
    df_motifs_bin['score_bin'] = 'Inactive'
    for score_lower_bound, score_upper_bound in zip(breaks[:-1], breaks[1:]):
        score_bin = 'Active: ' + str(score_lower_bound) + ('+' if score_upper_bound == float('inf') else ('-' + str(score_upper_bound)))
        score_bins.append(score_bin)
        df_temp = df_active.loc[df_active[score_col].between(score_lower_bound, score_upper_bound, inclusive='left')].copy()
        df_temp['score_bin'] = score_bin
        df_temp = df_temp.sample(sample_size) if len(df_temp)>1000 else df_temp
        df_motifs_bin = pd.concat([df_motifs_bin, df_temp], ignore_index=True)

    df_motifs_ct = df_motifs_bin.groupby('score_bin').size().reset_index()
    df_motifs_ct['score_bin'] = pd.Categorical(df_motifs_ct['score_bin'], categories=score_bins)
    df_motifs_ct.sort_values('score_bin', inplace=True)
    df_motifs_ct.rename(columns={0: 'motif_ct'}, inplace=True)
    print(df_motifs_ct)

    df_motifs_bin['score_bin'] = pd.Categorical(df_motifs_bin['score_bin'], categories=score_bins)
    return(df_motifs_bin)

def df_motifs_bin_2_df_read_type_frac_by_model(df_motifs_bin, bam, model, read_type_list = ['inactive', 'active'], model_type='1d'):
    valid_model_types = ('1d', '2d', 'oneHot')
    df_read_type_frac = pd.DataFrame()
    if model_type in valid_model_types:
        score_bins = df_motifs_bin['score_bin'].cat.categories
        for score_bin in score_bins:
            df_temp = df_motifs_bin.loc[df_motifs_bin['score_bin'] == score_bin]
            for index, row in tqdm(df_temp.iterrows(), total=len(df_temp)):
                # use model to predict read type
                if model_type == 'oneHot':
                    oneHotMat = np.expand_dims(motif_row_2_read_pairs_2_oneHotMat(bam, row), axis=3)
                    pred_probs = model.predict(oneHotMat, verbose=0)
                else:
                    mat_edits = motif_row_2_read_pairs_2_mat_edits(bam, row)
                    if model_type == '1d':
                        pred_probs = model.predict(mat_edits+1, verbose=0)
                    elif model_type == '2d':
                        mat_edits_2d = np.expand_dims(mat_edits_1d_to_2d(mat_edits), axis=3)
                        pred_probs = model.predict(mat_edits_2d, verbose=0)
                
                prob_sum = np.sum(pred_probs, axis=0)
                prob_frac = prob_sum/np.sum(prob_sum)
                pred_classes = np.argmax(pred_probs, axis=-1)
                read_ct = [sum(pred_classes == i) for i in range(len(read_type_list))]
                df_temp2 = pd.DataFrame({'score_bin': score_bin, 'chipseq_score': row['chipseq_score'], 
                                         'read_type': read_type_list, 
                                         'prob_sum': prob_sum, 'prob_frac': prob_frac/np.sum(prob_frac),
                                         'read_ct': read_ct, 'read_frac': [x/sum(read_ct) for x in read_ct]})

                df_read_type_frac = pd.concat([df_read_type_frac, df_temp2], ignore_index=True)
                    
        df_read_type_frac['score_bin'] = pd.Categorical(df_read_type_frac['score_bin'], categories=score_bins)
    else:
        print('Invalid model type! Must be 1d, 2d or oneHot.')
    return(df_read_type_frac)

def read2oneHotMat(read, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A']):
    edit_types = [x.split('2') for x in edit_types]
    flank_start = TF_center - flank_len
    flank_end = TF_center + flank_len

    # for a given read, one-hot encode the 4 bases and the 2 edit types
    # output as a list of lists for the flanking region
    aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
    q_seq = read.query_sequence.upper()
    oneHotMat = [[0] * (flank_len*2+1) for _ in range(6)]
    oneHot_dict = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'C2T': 4, 'G2A': 5}

    for q_pos, r_pos, ref in aligned_pairs:
        ref = ref.upper()
        if motif_strand == '+':
            i = r_pos - flank_start
        else:
            i = flank_end - r_pos

        if (i >= 0) and (i <= flank_len*2):
            oneHotMat[oneHot_dict[ref]][i] = 1
            alt = q_seq[q_pos]
            
            for ref2,alt2 in edit_types:
                if ref == ref2 and alt == alt2:
                    edit_type = ref2 + '2' + alt2
                    oneHotMat[oneHot_dict[edit_type]][i] = 1
                    break
    return(oneHotMat)

def merge_oneHotMats(oneHotMat_a, oneHotMat_b):
    merged_oneHotMat = [[]] * len(oneHotMat_a)
    for i in range(len(oneHotMat_a)):
        list_a = oneHotMat_a[i]
        list_b = oneHotMat_b[i]
        merged_list = merge_lists(list_a, list_b)
        merged_oneHotMat[i] = merged_list
    return(merged_oneHotMat)

def merge_read_pair_oneHotMats(read1, read2, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A']):
    if type(read1) != pysam.AlignedSegment:
        return read2oneHotMat(read2, TF_center, flank_len, motif_strand, edit_types)
    elif type(read2) != pysam.AlignedSegment:
        return read2oneHotMat(read1, TF_center, flank_len, motif_strand, edit_types)
    else:
        oneHotMat1 = read2oneHotMat(read1, TF_center, flank_len, motif_strand, edit_types)
        oneHotMat2 = read2oneHotMat(read2, TF_center, flank_len, motif_strand, edit_types)
        return(merge_oneHotMats(oneHotMat1, oneHotMat2))

def read_2_edit_channel(read, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A']):
    flank_start = TF_center - flank_len
    flank_end = TF_center + flank_len

    # for a given read, one-hot encode the 4 bases and the 2 edit types
    # output as a list of lists for the flanking region
    aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=True)
    q_seq = read.query_sequence.upper()
    edit_channel = [[0] * (flank_len*2+1) for _ in range(2)]
    oneHot_dict = {'C2T': 0, 'G2A': 1}
    rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    for q_pos, r_pos, ref in aligned_pairs:
        ref = ref.upper()
        alt = q_seq[q_pos]

        if (ref in rc_dict.keys()) and (alt in rc_dict.keys()):
            if motif_strand == '+':
                i = r_pos - flank_start
            else:
                i = flank_end - r_pos
                ref = rc_dict[ref]
                alt = rc_dict[alt]

            if (i >= 0) and (i <= flank_len*2):
                edit_type = ref+'2'+alt
                if edit_type in edit_types:
                    edit_channel[oneHot_dict[edit_type]][i] = 1
    return(edit_channel)

def merge_read_pair_edit_channel(read1, read2, TF_center, flank_len, motif_strand, edit_types = ['C2T', 'G2A']):
    if type(read1) != pysam.AlignedSegment:
        return read_2_edit_channel(read2, TF_center, flank_len, motif_strand, edit_types)
    elif type(read2) != pysam.AlignedSegment:
        return read_2_edit_channel(read1, TF_center, flank_len, motif_strand, edit_types)
    else:
        edit_channel1 = read_2_edit_channel(read1, TF_center, flank_len, motif_strand, edit_types)
        edit_channel2 = read_2_edit_channel(read2, TF_center, flank_len, motif_strand, edit_types)
        return(merge_oneHotMats(edit_channel1, edit_channel2))

def motif_row_2_base_channel(fasta, row, flank_len=150):
    ref_seq = fasta.fetch(row['chr'], row['center']-flank_len, row['center']+flank_len+1)
    if row['strand'] == '-':
        ref_seq = Seq(ref_seq)
        ref_seq = ref_seq.reverse_complement()

    base_channel = [[0] * (flank_len*2+1) for _ in range(4)]
    oneHot_dict = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    for i in range(len(ref_seq)):
        ref = ref_seq[i]
        if ref in oneHot_dict.keys():
            base_channel[oneHot_dict[ref]][i] = 1
    return(base_channel)
    
def df_motifs_2_read_pairs_2_oneHotMat(bam, fasta, df_motifs, flank_len=150,
                                       center_trough_pos_list = [x for x in range(-10,11)],
                                       left_peak_pos_list = [x for x in range(-75,-9)],
                                       right_peak_pos_list = [x for x in range(10,76)],
                                       center_trough_covs_frac_thres = 0.5,
                                       peak_covs_frac_thres = 0.5,
                                       read_pair_thres = 0,
                                       return_score_col = None):
    oneHotMat = []
    scoreList = []
    for index, row in tqdm(df_motifs.iterrows(), total=len(df_motifs)):
        TF_center = row['center']
        strand = row['strand']
        df_read_pairs = motif_row_2_df_read_pairs(bam, row)

        if len(df_read_pairs):
            df_read_pairs['chr'] = row['chr']
            df_read_pairs['center'] = TF_center
            df_read_pairs['strand'] = strand

            # get overall coverage list, filter reads by coverage in feature regions
            df_read_pairs['covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
            df_read_pairs['center_trough_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
            df_read_pairs['left_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, left_peak_pos_list, flank_len))
            df_read_pairs['right_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, right_peak_pos_list, flank_len))
            df_read_pairs = df_read_pairs.loc[(df_read_pairs['center_trough_covs_frac']>=center_trough_covs_frac_thres) & 
                                            ((df_read_pairs['left_peak_covs_frac']>=peak_covs_frac_thres) | (df_read_pairs['right_peak_covs_frac']>=peak_covs_frac_thres))]

            if len(df_read_pairs) > read_pair_thres:
                # score list
                if return_score_col in row.index:
                    df_read_pairs['return_score'] = row[return_score_col]
                    scoreList += df_read_pairs['return_score'].to_list()

                # base channel
                base_channel = motif_row_2_base_channel(fasta, row, flank_len)
                
                # oneHotMat_all = base_channel + covs_channel + edit_channel
                df_read_pairs['oneHotMat_all'] = df_read_pairs.apply(lambda row2: base_channel + [row2['covs']] + merge_read_pair_edit_channel(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
                oneHotMat += df_read_pairs['oneHotMat_all'].tolist()

    if return_score_col is not None:
        return(np.asarray(oneHotMat), scoreList)
    else:
        return(np.asarray(oneHotMat))
    
def motif_row_2_read_pairs_2_oneHotMat(bam, fasta, row, flank_len=150,
                                      center_trough_pos_list = [x for x in range(-10,11)],
                                      left_peak_pos_list = [x for x in range(-75,-9)],
                                      right_peak_pos_list = [x for x in range(10,76)],
                                      center_trough_covs_frac_thres = 0.5,
                                      peak_covs_frac_thres = 0.5,
                                      read_pair_thres = 0,
                                      return_score_col = None):
    oneHotMat = []
    scoreList = []
    TF_center = row['center']
    strand = row['strand']
    df_read_pairs = motif_row_2_df_read_pairs(bam, row)

    if len(df_read_pairs):
        df_read_pairs['chr'] = row['chr']
        df_read_pairs['center'] = TF_center
        df_read_pairs['strand'] = strand

        # get overall coverage list, filter reads by coverage in feature regions
        df_read_pairs['covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
        df_read_pairs['center_trough_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
        df_read_pairs['left_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, left_peak_pos_list, flank_len))
        df_read_pairs['right_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, right_peak_pos_list, flank_len))
        df_read_pairs = df_read_pairs.loc[(df_read_pairs['center_trough_covs_frac']>=center_trough_covs_frac_thres) & 
                                        ((df_read_pairs['left_peak_covs_frac']>=peak_covs_frac_thres) | (df_read_pairs['right_peak_covs_frac']>=peak_covs_frac_thres))]

        if len(df_read_pairs) > read_pair_thres:
            # score list
            if return_score_col in row.index:
                df_read_pairs['return_score'] = row[return_score_col]
                scoreList += df_read_pairs['return_score'].to_list()

            # base channel
            base_channel = motif_row_2_base_channel(fasta, row, flank_len)
            
            # oneHotMat_all = base_channel + covs_channel + edit_channel
            df_read_pairs['oneHotMat_all'] = df_read_pairs.apply(lambda row2: base_channel + [row2['covs']] + merge_read_pair_edit_channel(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
            oneHotMat += df_read_pairs['oneHotMat_all'].tolist()

    if return_score_col is not None:
        return(np.asarray(oneHotMat), scoreList)
    else:
        return(np.asarray(oneHotMat))
    
def df_read_pairs_2_oneHotMat(fasta, row, df_read_pairs, flank_len=150):
    TF_center = row['center']
    strand = row['strand']
    oneHotMat = []

    # base channel
    base_channel = motif_row_2_base_channel(fasta, row, flank_len)
    
    # oneHotMat_all = base_channel + covs_channel + edit_channel
    df_read_pairs['oneHotMat_all'] = df_read_pairs.apply(lambda row2: base_channel + [row2['covs']] + merge_read_pair_edit_channel(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
    oneHotMat += df_read_pairs['oneHotMat_all'].tolist()

    return(oneHotMat)

def df_motif_2_read_pairs_2_oneHotMat_dict(bam, fasta, df_motifs, flank_len=150,
                                           center_trough_pos_list = [x for x in range(-10,11)],
                                           left_peak_pos_list = [x for x in range(-75,-9)],
                                           right_peak_pos_list = [x for x in range(10,76)],
                                           center_trough_covs_frac_thres = 0.5,
                                           peak_covs_frac_thres = 0.5,
                                           read_ct_thres = 400000):
    peak_pos_list = left_peak_pos_list + right_peak_pos_list
    edit_thres_dict = {'recently bound': {'trough': (0.2, 1,'both'), 'peak': (0.2, 1, 'both')},
                       'bound': {'trough': (0, 0.2,'left'), 'peak': (0.2, 1, 'both')},
                       'unbound': {'trough': (0, 1,'both'), 'peak': (0, 0.2, 'left')}}
    edit_type_list = ['C2T', 'G2A']
    editable_base_list = ['C', 'G']
    oneHotMat_dict = {'recently bound': [], 'bound': [], 'unbound': []}

    motif_ct = 0
    for index, row in tqdm(df_motifs.iterrows(), total=len(df_motifs)):
        TF_center = row['center']
        strand = row['strand']
        df_read_pairs = motif_row_2_df_read_pairs(bam, row)

        if len(df_read_pairs):
            df_read_pairs['chr'] = row['chr']
            df_read_pairs['center'] = TF_center
            df_read_pairs['strand'] = strand

            # get overall coverage list, filter reads by coverage in feature regions
            df_read_pairs['covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand), axis=1)
            df_read_pairs['center_trough_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
            df_read_pairs['left_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, left_peak_pos_list, flank_len))
            df_read_pairs['right_peak_covs_frac'] = df_read_pairs['covs'].apply(lambda x: covs_frac_in_range(x, right_peak_pos_list, flank_len))
            df_read_pairs = df_read_pairs.loc[(df_read_pairs['center_trough_covs_frac']>=center_trough_covs_frac_thres) & 
                                            ((df_read_pairs['left_peak_covs_frac']>=peak_covs_frac_thres) | (df_read_pairs['right_peak_covs_frac']>=peak_covs_frac_thres))]

            if len(df_read_pairs):
                for edit_type, editable_base in zip(edit_type_list, editable_base_list):
                    # get edits list, editable coverage list, and calculate feature region edit fraction
                    df_read_pairs['edits'] = df_read_pairs.apply(lambda row2: merge_read_pair_edits(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand, edit_types=[edit_type]), axis=1)
                    df_read_pairs['editable_covs'] = df_read_pairs.apply(lambda row2: merge_read_pair_covs(row2['r1_reads'], row2['r2_reads'], TF_center, flank_len, strand, editable_pos_only=True, editable_bases=[editable_base]), axis=1)
                    df_read_pairs['trough_edits_frac'] = df_read_pairs.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
                    df_read_pairs['peak_edits_frac'] = df_read_pairs.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
                    
                    for read_type in edit_thres_dict:
                        if len(oneHotMat_dict[read_type]) < read_ct_thres:
                            peak_edit_thres1, peak_edit_thres2, peak_edit_inclusive = edit_thres_dict[read_type]['peak']
                            trough_edit_thres1, trough_edit_thres2, trough_edit_inclusive = edit_thres_dict[read_type]['trough']
                            df_temp = df_read_pairs.loc[(df_read_pairs['trough_edits_frac'].between(trough_edit_thres1, trough_edit_thres2, inclusive=trough_edit_inclusive)) & 
                                                        (df_read_pairs['peak_edits_frac'].between(peak_edit_thres1, peak_edit_thres2, inclusive=peak_edit_inclusive))].copy()
                            if len(df_temp):
                                oneHotMat_dict[read_type] += df_read_pairs_2_oneHotMat(fasta, row, df_temp, flank_len=flank_len)

            if len(oneHotMat_dict['recently bound']) >= read_ct_thres and len(oneHotMat_dict['bound']) >= read_ct_thres and len(oneHotMat_dict['unbound']) >= read_ct_thres:
                break
        
        motif_ct += 1
        if motif_ct % 1000 == 0:
            print(f"Unbound reads: {len(oneHotMat_dict['unbound'])}, bound reads: {len(oneHotMat_dict['bound'])}, recently bound reads: {len(oneHotMat_dict['recently bound'])}")

    for edit_type in oneHotMat_dict:
        oneHotMat_dict[edit_type] = np.asarray(oneHotMat_dict[edit_type])
        print(f'{edit_type}: {oneHotMat_dict[edit_type].shape}')
    
    return(oneHotMat_dict)


def oneHotMat_2_edit_frac(oneHotMat):
    flank_len = int((oneHotMat.shape[2] - 1)/2)
    edits = np.zeros(shape=oneHotMat.shape[2])
    covs = np.zeros(shape=oneHotMat.shape[2])

    for i in range(oneHotMat.shape[0]):
        C2T_edits = oneHotMat[i,5,:]
        G2A_edits = oneHotMat[i,6,:]
        covs_G = ((oneHotMat[i,2,:] == 1) & (oneHotMat[i,4,:] == 1)).astype(int)
        covs_C = ((oneHotMat[i,3,:] == 1) & (oneHotMat[i,4,:] == 1)).astype(int)

        read_edit_type = get_edit_type(sum(C2T_edits), sum(G2A_edits))
        if read_edit_type == 'C2T_dominant':
            edits += C2T_edits
            covs += covs_C
        elif read_edit_type == 'G2A_dominant':
            edits += G2A_edits
            covs += covs_G
        else:
            edits += (C2T_edits + G2A_edits)
            covs += (covs_C + covs_G)
            
    edit_frac = edits / covs
    df_edit_frac = pd.DataFrame({'relative_pos': [x for x in range(-flank_len, flank_len+1)], 'edit_frac': edit_frac})
    return(df_edit_frac)

def oneHotMat_2_edit_frac_by_read_type_labels(oneHotMat, read_type_labels, read_type_list = ['unbound', 'bound', 'recently bound']):
    df_composite_edit_frac = pd.DataFrame()
    for i in range(len(read_type_list)):
        oneHotMat_label = oneHotMat[read_type_labels == i,:,:]
        df_temp = oneHotMat_2_edit_frac(oneHotMat_label)
        df_temp['read_type'] = read_type_list[i]
        df_composite_edit_frac = pd.concat([df_composite_edit_frac, df_temp])
    return(df_composite_edit_frac)
        

def create_model(input_dim, conv_layer_ct=3, dense_layer_ct=2, 
                 conv_starting_node_ct=16, filter_size=12, pool_size=(1,2), 
                 dense_starting_node_ct=32, 
                 output_node_ct=2,
                 dropout_rate=0.25, l2_reg_rate=0.01, learn_rate=0.0001):
    l2_reg = tf.keras.regularizers.l2(l2_reg_rate)

    model = tf.keras.Sequential()
    model.add(tf.keras.Input(shape=input_dim))

    for i in range(conv_layer_ct):
        if i == 0:
            conv_filter_size = (input_dim[0],filter_size)
        else:
            conv_filter_size = (1,filter_size)
        
        model.add(tf.keras.layers.Conv2D(conv_starting_node_ct*(i+1), conv_filter_size, kernel_regularizer=l2_reg))
        model.add(tf.keras.layers.BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(tf.keras.layers.MaxPool2D(pool_size))
        model.add(tf.keras.layers.Dropout(dropout_rate))

    model.add(tf.keras.layers.Flatten())

    for i in range(dense_layer_ct):
        model.add(tf.keras.layers.Dense(int(dense_starting_node_ct/(2**i)), kernel_regularizer=l2_reg))
        model.add(tf.keras.layers.BatchNormalization())
        model.add(tf.keras.layers.Activation('relu'))
        model.add(tf.keras.layers.Dropout(dropout_rate))

    model.add(tf.keras.layers.Dense(output_node_ct, activation='softmax'))

    # Compile the model
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=learn_rate), 
                loss='categorical_crossentropy', metrics=['accuracy'])
    return(model)

def get_model_history(history):
    df_history = pd.DataFrame(history.history)
    df_history['epoch'] = [x for x in range(len(df_history))]
    df_history.rename(columns={'accuracy': 'accuracy_train', 'loss': 'loss_train', 
                            'val_accuracy': 'accuracy_validation', 'val_loss': 'loss_validation'}, inplace=True)
    df_history = pd.melt(df_history, id_vars=['epoch'], value_vars=['accuracy_train', 'loss_train', 'accuracy_validation', 'loss_validation'])
    df_history['type'] = df_history['variable'].apply(lambda x: x.split('_')[0])
    df_history['val_type'] = df_history['variable'].apply(lambda x: x.split('_')[1])
    return(df_history)


def df_reads_2_composite_edit_frac(df_reads, flank_len, edit_types = ['C2T', 'G2A']):
    editable_bases = [x.split('2')[0] for x in edit_types]
    edit_ct = np.sum([0] * (flank_len*2+1))
    total_ct = np.sum([0] * (flank_len*2+1))
    
    # for C2T and G2A dominant reads, edit fractions are calculated using only the dominant edit types
    for edit_type, editable_base in zip(edit_types, editable_bases):
        df_temp = df_reads.loc[df_reads['read_edit_type'] == (edit_type + '_dominant')].copy()
        if len(df_temp):
            if (edit_type + '_edits') in df_temp.columns:
                edit_mat = np.array(df_temp[edit_type + '_edits'].tolist())
            else:
                edit_mat = np.array(df_temp.apply(lambda row2: read2edits(row2['reads'], row2['center'], flank_len, row2['strand'], edit_types=[edit_type]), axis=1).tolist())
            editable_covs_mat = np.array(df_temp.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, row2['strand'], editable_pos_only=True, editable_bases=[editable_base]), axis=1).tolist())

            edit_ct += np.sum(edit_mat, axis=0)
            total_ct += np.sum(editable_covs_mat, axis=0)
    
    # for mixed_edits and no_edits reads, edit fractions are calculated using both edit types
    df_temp = df_reads.loc[df_reads['read_edit_type'].isin(['mixed_edits', 'no_edits'])].copy()
    if len(df_temp):
        for edit_type, editable_base in zip(edit_types, editable_bases):
            if (edit_type + '_edits') in df_temp.columns:
                edit_mat = np.array(df_temp[edit_type + '_edits'].tolist())
            else:
                edit_mat = np.array(df_temp.apply(lambda row2: read2edits(row2['reads'], row2['center'], flank_len, row2['strand'], edit_types=[edit_type]), axis=1).tolist())
            editable_covs_mat = np.array(df_temp.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, row2['strand'], editable_pos_only=True, editable_bases=[editable_base]), axis=1).tolist())
            edit_ct += np.sum(edit_mat, axis=0)
            total_ct += np.sum(editable_covs_mat, axis=0)
    
    df_baseCt = pd.DataFrame({'relative_pos': [i for i in range(-flank_len, flank_len+1)],
                              'edit_ct': edit_ct,
                              'total_ct': total_ct})
    df_baseCt['edit_frac'] = df_baseCt['edit_ct'] / df_baseCt['total_ct']
    return(df_baseCt)


def motif_row_2_df_reads(bam, row):
    chr = row['chr']
    flank_start = row['flank_start']
    flank_end = row['flank_end']

    r_names = []
    r_reads = []
    for read in bam.fetch(chr, flank_start, flank_end):
        read_name = read.query_name
        r_names.append(read_name)
        r_reads.append(read) 
    
    df_reads = pd.DataFrame({'names': r_names, 'reads': r_reads})
    return(df_reads)

def get_edit_type(C2T_edits, G2A_edits, edit_per_read_thres=1, edit_majority_thres=0.6):
    # assign read label according to dominant edits
    total_edits = C2T_edits + G2A_edits
    edit_type = 'no_edits'
    if total_edits > edit_per_read_thres:
        if C2T_edits/total_edits > edit_majority_thres:
            edit_type = 'C2T_dominant'
        elif G2A_edits/total_edits > edit_majority_thres:
            edit_type = 'G2A_dominant'
        else:
            edit_type = 'mixed_edits'
    return edit_type


def df_motifs_2_df_reads_2_composite_edit_frac(bam, df_motifs, flank_len=300):    
    ct=0
    df_baseCt_grouped = pd.DataFrame()
    for index, row in df_motifs.iterrows():
        TF_center = row['center']
        strand = row['strand']
        df_reads = motif_row_2_df_reads(bam, row)

        if len(df_reads):
            df_reads['chr'] = row['chr']
            df_reads['center'] = TF_center
            df_reads['strand'] = strand

            # get edits list
            df_reads['C2T_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['C2T']), axis=1)
            df_reads['G2A_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['G2A']), axis=1)
            df_reads['read_edit_type'] = df_reads.apply(lambda row2: get_edit_type(sum(row2['C2T_edits']), sum(row2['G2A_edits'])), axis=1)
        
            df_baseCt = df_reads_2_composite_edit_frac(df_reads, flank_len)
            df_baseCt_grouped = pd.concat([df_baseCt_grouped, df_baseCt[['relative_pos', 'edit_ct', 'total_ct']]], axis=0, ignore_index=True)

        ct += 1
        if ct % 1000 == 0:
            df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'}).reset_index()
            print(f'Motifs processed: {ct}')

    df_baseCt_grouped = df_baseCt_grouped.groupby(['relative_pos']).agg({'edit_ct': 'sum', 'total_ct': 'sum'})
    df_baseCt_grouped = df_baseCt_grouped.reset_index()
    df_baseCt_grouped['edit_frac'] = df_baseCt_grouped['edit_ct'] / df_baseCt_grouped['total_ct']
    return(df_baseCt_grouped)


def prepare_df_reads_3readType(df_reads, row, flank_len = 150, 
                               center_trough_pos_list = [x for x in range(-10, 11)],
                               center_trough_covs_frac_thres = 0.5):
        chrom = row['chr']
        TF_center = row['center']
        strand = row['strand']
        df_reads['chr'] = chrom
        df_reads['center'] = TF_center
        df_reads['strand'] = strand

        # get edits list
        df_reads['C2T_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['C2T']), axis=1)
        df_reads['G2A_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['G2A']), axis=1)
        df_reads['read_edit_type'] = df_reads.apply(lambda row2: get_edit_type(sum(row2['C2T_edits']), sum(row2['G2A_edits'])), axis=1)
        df_reads['edits'] = df_reads.apply(lambda row2: row2['C2T_edits'] if row2['read_edit_type']=='C2T_dominant' else row2['G2A_edits'] if row2['read_edit_type']=='G2A_dominant' else [x + y for x, y in zip(row2['C2T_edits'], row2['G2A_edits'])], axis=1)
        
        # get coverage list
        df_reads['covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand), axis=1)
        df_reads['editable_C_covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand, editable_pos_only=True, editable_bases=['C']), axis=1)
        df_reads['editable_G_covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand, editable_pos_only=True, editable_bases=['G']), axis=1)
        df_reads['editable_covs'] = df_reads.apply(lambda row2: row2['editable_C_covs'] if row2['read_edit_type']=='C2T_dominant' else row2['editable_G_covs'] if row2['read_edit_type']=='G2A_dominant' else [x + y for x, y in zip(row2['editable_C_covs'], row2['editable_G_covs'])], axis=1)

        # filter reads by coverage in feature regions
        df_reads['center_trough_covs_frac'] = df_reads['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
        df_reads = df_reads.loc[(df_reads['center_trough_covs_frac']>=center_trough_covs_frac_thres)]

        return(df_reads)


def df_reads_2_oneHotMat(fasta, row, df_reads, flank_len=150, np_array=True):
    TF_center = row['center']
    strand = row['strand']
    oneHotMat = []

    # base channel
    base_channel = motif_row_2_base_channel(fasta, row, flank_len)
    
    # oneHotMat_all = base_channel + covs_channel + edit_channel
    df_reads['oneHotMat_all'] = df_reads.apply(lambda row2: base_channel + [row2['covs']] + read_2_edit_channel(row2['reads'], TF_center, flank_len, strand), axis=1)
    oneHotMat += df_reads['oneHotMat_all'].tolist()

    if np_array:
        return(np.array(oneHotMat))
    else:
        return(oneHotMat)


def df_motif_2_df_reads_2_oneHotMat_dict(bam, fasta, df_motifs, flank_len=150,
                                         center_trough = [-10,10], left_peak = [-75, -10], right_peak = [10, 75],
                                         center_trough_covs_frac_thres = 0.5, trough_edit_thres = 0.2, peak_edit_thres = 0.2, global_edit_thres = 0.2,
                                         read_ct_thres = 400000):
    
    center_trough_pos_list = [x for x in range(center_trough[0], center_trough[1]+1)]
    left_peak_pos_list = [x for x in range(left_peak[0], left_peak[1]+1)]
    right_peak_pos_list = [x for x in range(right_peak[0], right_peak[1]+1)]
    peak_pos_list = left_peak_pos_list + right_peak_pos_list

    edit_thres_dict = {'recently bound': {'trough': [trough_edit_thres, 1,'both'], 'peak': [peak_edit_thres, 1, 'both'], 'global': [0, 1, 'both']},
                       'bound': {'trough': [0, trough_edit_thres,'left'], 'peak': [peak_edit_thres, 1, 'both'], 'global': [0, 1, 'both']},
                       'unbound': {'trough': [0, 1,'both'], 'peak': [0, peak_edit_thres, 'both'], 'global': [0, global_edit_thres, 'left']}}
    oneHotMat_dict = {'recently bound': [], 'bound': [], 'unbound': []}

    motif_ct = 0
    for index, row in df_motifs.sample(frac=1).iterrows():
        df_reads = motif_row_2_df_reads(bam, row)

        if len(df_reads):
            df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

            if len(df_reads):
                df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
                df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
                df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)

                for read_type in edit_thres_dict:
                    if len(oneHotMat_dict[read_type]) <= read_ct_thres:
                        peak_edit_thres1, peak_edit_thres2, peak_edit_inclusive = edit_thres_dict[read_type]['peak']
                        trough_edit_thres1, trough_edit_thres2, trough_edit_inclusive = edit_thres_dict[read_type]['trough']
                        global_edit_thres1, global_edit_thres2, global_edit_inclusive = edit_thres_dict[read_type]['global']
                        df_temp = df_reads.loc[(df_reads['trough_edits_frac'].between(trough_edit_thres1, trough_edit_thres2, inclusive=trough_edit_inclusive)) & 
                                                    (df_reads['peak_edits_frac'].between(peak_edit_thres1, peak_edit_thres2, inclusive=peak_edit_inclusive)) &
                                                    (df_reads['global_edits_frac'].between(global_edit_thres1, global_edit_thres2, inclusive=global_edit_inclusive))].copy()
                        if len(df_temp):
                            oneHotMat_dict[read_type] += df_reads_2_oneHotMat(fasta, row, df_temp, flank_len=flank_len, np_array=False)

        if (len(oneHotMat_dict['recently bound']) >= read_ct_thres) and (len(oneHotMat_dict['bound']) >= read_ct_thres) and (len(oneHotMat_dict['unbound']) >= read_ct_thres):
            print("All labels exceed read count threshold.")
            break
    
        motif_ct += 1
        if motif_ct % 1000 == 0:
            print(f"Unbound reads: {len(oneHotMat_dict['unbound'])}, bound reads: {len(oneHotMat_dict['bound'])}, recently bound reads: {len(oneHotMat_dict['recently bound'])}")

    for edit_type in oneHotMat_dict:
        oneHotMat_dict[edit_type] = np.asarray(oneHotMat_dict[edit_type])
        print(f'{edit_type}: {oneHotMat_dict[edit_type].shape}')
    
    return(oneHotMat_dict)


def motif_row_2_df_reads_2_oneHotMat_dict(bam, fasta, row, flank_len=150,
                                          center_trough = [-10,10], left_peak = [-75, -10], right_peak = [10, 75],
                                          center_trough_covs_frac_thres = 0.5, trough_edit_thres = 0.2, peak_edit_thres = 0.2, global_edit_thres = 0.2):
    
    center_trough_pos_list = [x for x in range(center_trough[0], center_trough[1]+1)]
    left_peak_pos_list = [x for x in range(left_peak[0], left_peak[1]+1)]
    right_peak_pos_list = [x for x in range(right_peak[0], right_peak[1]+1)]
    peak_pos_list = left_peak_pos_list + right_peak_pos_list

    edit_thres_dict = {'recently bound': {'trough': [trough_edit_thres, 1,'both'], 'peak': [peak_edit_thres, 1, 'both'], 'global': [0, 1, 'both']},
                       'bound': {'trough': [0, trough_edit_thres,'left'], 'peak': [peak_edit_thres, 1, 'both'], 'global': [0, 1, 'both']},
                       'unbound': {'trough': [0, 1,'both'], 'peak': [0, peak_edit_thres, 'both'], 'global': [0, global_edit_thres, 'left']}}
    oneHotMat_dict = {'recently bound': [], 'bound': [], 'unbound': []}

    df_reads = motif_row_2_df_reads(bam, row)

    if len(df_reads):
        df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

        if len(df_reads):
            df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
            df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
            df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)

            for read_type in edit_thres_dict:
                peak_edit_thres1, peak_edit_thres2, peak_edit_inclusive = edit_thres_dict[read_type]['peak']
                trough_edit_thres1, trough_edit_thres2, trough_edit_inclusive = edit_thres_dict[read_type]['trough']
                global_edit_thres1, global_edit_thres2, global_edit_inclusive = edit_thres_dict[read_type]['global']
                df_temp = df_reads.loc[(df_reads['trough_edits_frac'].between(trough_edit_thres1, trough_edit_thres2, inclusive=trough_edit_inclusive)) & 
                                            (df_reads['peak_edits_frac'].between(peak_edit_thres1, peak_edit_thres2, inclusive=peak_edit_inclusive)) &
                                            (df_reads['global_edits_frac'].between(global_edit_thres1, global_edit_thres2, inclusive=global_edit_inclusive))].copy()
                if len(df_temp):
                    oneHotMat_dict[read_type] += df_reads_2_oneHotMat(fasta, row, df_temp, flank_len=flank_len, np_array=False)

    for edit_type in oneHotMat_dict:
        oneHotMat_dict[edit_type] = np.asarray(oneHotMat_dict[edit_type])
    
    return(oneHotMat_dict)


def spearmanr_wrapper(list_a, list_b, prefix=''):
    cor_val, p_val = spearmanr(list_a, list_b, nan_policy='omit')
    return {f'{prefix}rho': cor_val, f'{prefix}pval': p_val}

def pearsonr_wrapper(list_a, list_b, prefix=''):
    array_a = np.array(list_a)
    array_b = np.array(list_b)
    
    # Create a mask to filter out NaN values in both arrays
    mask = ~np.isnan(array_a) & ~np.isnan(array_b)
    clean_a = array_a[mask]
    clean_b = array_b[mask]

    cor_val, p_val = pearsonr(clean_a, clean_b)
    return {f'{prefix}r': cor_val, f'{prefix}pval': p_val}


def df_reads_2_predicted_labels(df_reads, model, fasta, row, flank_len=150, read_type_list = ['unbound', 'bound', 'recently bound']):
    # generate model input matrix
    oneHotMat = df_reads_2_oneHotMat(fasta, row, df_reads, flank_len)

    # use model to predict read type
    oneHotMat_conv = np.expand_dims(oneHotMat, axis=3)
    pred_probs = model.predict(oneHotMat_conv, verbose=0)
    pred_classes = np.argmax(pred_probs, axis=1)

    return [read_type_list[i] for i in pred_classes]
    
def df_reads_2_predicted_probs(df_reads, model, fasta, row, flank_len=150, read_type_list = ['unbound', 'bound', 'recently bound']):
    # generate model input matrix
    oneHotMat = df_reads_2_oneHotMat(fasta, row, df_reads, flank_len)

    # use model to predict read type
    oneHotMat_conv = np.expand_dims(oneHotMat, axis=3)
    pred_probs = model.predict(oneHotMat_conv, verbose=0)
    pred_classes = np.argmax(pred_probs, axis=1)

    df_pred_probs = pd.DataFrame(pred_probs, columns=['read_prob ' + read_type for read_type in read_type_list])
    df_pred_probs['read_type'] = [read_type_list[i] for i in pred_classes]
    df_pred_probs['names'] = df_reads['names'].tolist()
    
    return df_pred_probs

def filter_df_motifs(df_motifs, fimo_pVal_thres=-np.log10(0.0001), chipseq_qVal_thres=-np.log10(0.05), read_ct_thres=50):
    df_active_pool = df_motifs.loc[(df_motifs['fimo_pVal'] >= fimo_pVal_thres) & ((df_motifs['motif_type'] == 'Active') & (df_motifs['chipseq_qVal'] > chipseq_qVal_thres)) & (df_motifs['motif_read_ct'] >= read_ct_thres)].copy()
    df_inactive_pool = df_motifs.loc[(df_motifs['fimo_pVal'] >= fimo_pVal_thres) & ((df_motifs['motif_type'] == 'Inactive') & (df_motifs['chipseq_qVal'] < 0)) & (df_motifs['motif_read_ct'] >= read_ct_thres)].copy()
    df_motifs_filtered = pd.concat([df_active_pool, df_inactive_pool], ignore_index=True)
    return df_motifs_filtered

def adjust_motif_center(row, df_motif_features, cell_type):
    row_TF = df_motif_features.loc[(df_motif_features['TF'] == row['TF']) & (df_motif_features['cell_type'] == cell_type)].iloc[0]
    center_adj_TF = row_TF[['l_peak_l', 'l_peak_r', 'r_peak_l', 'r_peak_r', 'center_l', 'center_r']].mean()
    center_TF = row[['start', 'end']].mean()
    center_TF = (center_TF + center_adj_TF) if (row['strand'] == '+') else (center_TF- center_adj_TF)
    return center_TF

def agg_by_bin_window(df_cobinding_bound, bin_col='delta_oe_prob', max_center_dist=100, one_side_window_size=1):
    center_dist_list = []
    observed_expected_delta_ma_list = []
    for center_dist in range(-max_center_dist, max_center_dist+1):
        center_dist_list.append(center_dist)
        observed_expected_delta_ma = df_cobinding_bound.loc[df_cobinding_bound['center_dist'].between(center_dist-one_side_window_size-0.5, center_dist+one_side_window_size+0.5, inclusive='left'), bin_col].median()
        observed_expected_delta_ma_list.append(observed_expected_delta_ma)
    df_cobinding_bound_ma = pd.DataFrame({'center_dist': center_dist_list, (bin_col + '_by_pos'): observed_expected_delta_ma_list})
    return df_cobinding_bound_ma


def calc_motif_features(df_composite_editFrac, flank_len=150):
    df_composite_editFrac['edit_frac_ma'] = df_composite_editFrac['edit_frac'].rolling(window=15, min_periods=1, center=True).mean()

    # calculate key edit fraction values
    min_center_editFrac_range = [-10, 10]
    max_l_peak_editFrac_ma_range = [-75, -14]
    max_r_peak_editFrac_ma_range = [14, 75]

    min_center_editFrac = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(min_center_editFrac_range[0], min_center_editFrac_range[1]), 'edit_frac'].min()
    max_l_peak_editFrac_ma = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(max_l_peak_editFrac_ma_range[0], max_l_peak_editFrac_ma_range[1]), 'edit_frac_ma'].max()
    max_r_peak_editFrac_ma = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(max_r_peak_editFrac_ma_range[0], max_r_peak_editFrac_ma_range[1]), 'edit_frac_ma'].max()
    min_peak_editFrac = min(max_l_peak_editFrac_ma, max_r_peak_editFrac_ma)
    mean_flankLimit_editFrac = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].isin([-flank_len, flank_len]), 'edit_frac_ma'].mean()
    # mean_flankLimit_editFrac = df_composite_editFrac['edit_frac_ma'].values.mean()

    # define footprint
    center_detect_range = [-12, 12]
    trough_detect_thres = 0.55

    df_composite_editFrac_temp = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(center_detect_range[0], center_detect_range[1])].copy()
    df_composite_editFrac_temp['ratio'] = (df_composite_editFrac_temp['edit_frac'] - min_center_editFrac) / (min_peak_editFrac - min_center_editFrac)
    center_pos_list = df_composite_editFrac_temp.loc[df_composite_editFrac_temp['ratio'] < trough_detect_thres, 'relative_pos'].tolist()

    # define left/right peak
    l_peak_detect_range = [-100, max(min(center_pos_list)-1, -100)]
    r_peak_detect_range = [min(max(center_pos_list)+1, 100), 100]
    peak_detect_thres = 2/3

    df_composite_editFrac_temp = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(l_peak_detect_range[0], l_peak_detect_range[1])].copy()
    df_composite_editFrac_temp['ratio'] = (df_composite_editFrac_temp['edit_frac'] - mean_flankLimit_editFrac) / (max_l_peak_editFrac_ma - mean_flankLimit_editFrac)
    l_peak_pos_list = df_composite_editFrac_temp.loc[df_composite_editFrac_temp['ratio'] > peak_detect_thres, 'relative_pos'].tolist()

    df_composite_editFrac_temp = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].between(r_peak_detect_range[0], r_peak_detect_range[1])].copy()
    df_composite_editFrac_temp['ratio'] = (df_composite_editFrac_temp['edit_frac'] - mean_flankLimit_editFrac) / (max_r_peak_editFrac_ma - mean_flankLimit_editFrac)
    r_peak_pos_list = df_composite_editFrac_temp.loc[df_composite_editFrac_temp['ratio'] > peak_detect_thres, 'relative_pos'].tolist()

    center_mean_editFrac = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].isin(center_pos_list), 'edit_frac'].mean()
    l_peak_mean_editFrac = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].isin(l_peak_pos_list), 'edit_frac'].mean()
    r_peak_mean_editFrac = df_composite_editFrac.loc[df_composite_editFrac['relative_pos'].isin(r_peak_pos_list), 'edit_frac'].mean()

    df_motif_features = pd.DataFrame({'relative_pos': center_pos_list + l_peak_pos_list + r_peak_pos_list,
                                    'feature': ['footprint'] * len(center_pos_list) + ['left_peak'] * len(l_peak_pos_list) + ['right_peak'] * len(r_peak_pos_list),
                                    'mean_editFrac': [center_mean_editFrac] * len(center_pos_list) + [l_peak_mean_editFrac] * len(l_peak_pos_list) + [r_peak_mean_editFrac] * len(r_peak_pos_list)})
    return df_motif_features

def calc_motif_features_v2(df_editFrac_active, df_editFrac_inactive, flank_len=150):
    df_editFrac_active['edit_frac_ma'] = df_editFrac_active['edit_frac'].rolling(window=15, min_periods=1, center=True).mean()

    # calculate key edit fraction values
    min_center_editFrac_range = [-10, 10]
    max_l_peak_editFrac_ma_range = [-75, -14]
    max_r_peak_editFrac_ma_range = [14, 75]

    min_center_editFrac = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(min_center_editFrac_range[0], min_center_editFrac_range[1]), 'edit_frac'].min()
    max_l_peak_editFrac_ma = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(max_l_peak_editFrac_ma_range[0], max_l_peak_editFrac_ma_range[1]), 'edit_frac_ma'].max()
    max_r_peak_editFrac_ma = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(max_r_peak_editFrac_ma_range[0], max_r_peak_editFrac_ma_range[1]), 'edit_frac_ma'].max()
    min_peak_editFrac = min(max_l_peak_editFrac_ma, max_r_peak_editFrac_ma)

    # define footprint
    center_detect_range = [-12, 12]
    trough_detect_thres = 0.55

    df_editFrac_active_temp = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(center_detect_range[0], center_detect_range[1])].copy()
    df_editFrac_active_temp['ratio'] = (df_editFrac_active_temp['edit_frac'] - min_center_editFrac) / (min_peak_editFrac - min_center_editFrac)
    center_pos_list = df_editFrac_active_temp.loc[df_editFrac_active_temp['ratio'] < trough_detect_thres, 'relative_pos'].tolist()

    # define left/right peak
    l_peak_detect_range = [-100, max(min(center_pos_list)-1, -100)]
    r_peak_detect_range = [min(max(center_pos_list)+1, 100), 100]
    peak_detect_thres = 2/3

    df_editFrac_active_temp = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(l_peak_detect_range[0], l_peak_detect_range[1])].copy()
    df_editFrac_active_temp['ratio'] = (df_editFrac_active_temp['edit_frac'] - mean_flankLimit_editFrac) / (max_l_peak_editFrac_ma - mean_flankLimit_editFrac)
    l_peak_pos_list = df_editFrac_active_temp.loc[df_editFrac_active_temp['ratio'] > peak_detect_thres, 'relative_pos'].tolist()

    df_editFrac_active_temp = df_editFrac_active.loc[df_editFrac_active['relative_pos'].between(r_peak_detect_range[0], r_peak_detect_range[1])].copy()
    df_editFrac_active_temp['ratio'] = (df_editFrac_active_temp['edit_frac'] - mean_flankLimit_editFrac) / (max_r_peak_editFrac_ma - mean_flankLimit_editFrac)
    r_peak_pos_list = df_editFrac_active_temp.loc[df_editFrac_active_temp['ratio'] > peak_detect_thres, 'relative_pos'].tolist()

    # calculate mean edit fraction for active motif features
    center_mean_editFrac = df_editFrac_active.loc[df_editFrac_active['relative_pos'].isin(center_pos_list), 'edit_frac'].mean()
    l_peak_mean_editFrac = df_editFrac_active.loc[df_editFrac_active['relative_pos'].isin(l_peak_pos_list), 'edit_frac'].mean()
    r_peak_mean_editFrac = df_editFrac_active.loc[df_editFrac_active['relative_pos'].isin(r_peak_pos_list), 'edit_frac'].mean()

    # calculate global mean edit fraction and std for inactive motifs
    inactive_motif_editFrac_mean = df_editFrac_inactive['edit_frac'].mean()
    inactive_motif_editFrac_std = df_editFrac_inactive['edit_frac'].std() 

    df_motif_features = pd.DataFrame({'center_pos_list': [center_pos_list],
                                      'l_peak_pos_list': [l_peak_pos_list],
                                      'r_peak_pos_list': [r_peak_pos_list],
                                      'footprint': [center_mean_editFrac],
                                      'left_peak': [l_peak_mean_editFrac],
                                      'right_peak': [r_peak_mean_editFrac],
                                      'inactive_mean': [inactive_motif_editFrac_mean],
                                      'inactive_std': [inactive_motif_editFrac_std]})
    return df_motif_features

def expand_active_motif_features(df_motif_features):
    df_active_motif_features_expanded_all = pd.DataFrame()
    for index, row in df_motif_features.iterrows():
        center_pos_list = row['center_pos_list']
        l_peak_pos_list = row['l_peak_pos_list']
        r_peak_pos_list = row['r_peak_pos_list']
        center_mean_editFrac = row['footprint']
        l_peak_mean_editFrac = row['left_peak']
        r_peak_mean_editFrac = row['right_peak']

        df_active_motif_features_expanded = pd.DataFrame({'relative_pos': center_pos_list + l_peak_pos_list + r_peak_pos_list,
                                                 'feature': ['footprint'] * len(center_pos_list) + ['left_peak'] * len(l_peak_pos_list) + ['right_peak'] * len(r_peak_pos_list),
                                                 'mean_editFrac': [center_mean_editFrac] * len(center_pos_list) + [l_peak_mean_editFrac] * len(l_peak_pos_list) + [r_peak_mean_editFrac] * len(r_peak_pos_list)})
        df_active_motif_features_expanded_all = pd.concat([df_active_motif_features_expanded_all, df_active_motif_features_expanded], ignore_index=True)

    return df_active_motif_features_expanded_all


def df_motif_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, df_motifs, 
                                            center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                            trough_thres, peak_thres, global_thres,
                                            flank_len=150, center_trough_covs_frac_thres = 0.5, 
                                            read_ct_thres = 400000):
    
    peak_pos_list = left_peak_pos_list + right_peak_pos_list
    oneHotMat_dict = {'recently bound': [], 'bound': [], 'unbound': []}

    motif_ct = 0
    for index, row in df_motifs.sample(frac=1).iterrows():
        df_reads = motif_row_2_df_reads(bam, row)

        if len(df_reads):
            df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

            if len(df_reads):
                df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
                df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
                df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)
                df_reads = df_reads.dropna(subset=['trough_edits_frac', 'peak_edits_frac', 'global_edits_frac'])

                # label unbound bound, and recently bound reads
                df_unbound_reads = df_reads.loc[(df_reads['peak_edits_frac'] < peak_thres) & (df_reads['global_edits_frac'] < global_thres)].copy()
                df_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (df_reads['trough_edits_frac'] < trough_thres)].copy()
                df_recently_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (~df_reads.index.isin(df_bound_reads.index))].copy()

                # calculate one hot matrix for each read type
                if len(oneHotMat_dict['unbound']) <= read_ct_thres and len(df_unbound_reads) > 0:
                    oneHotMat_dict['unbound'] += df_reads_2_oneHotMat(fasta, row, df_unbound_reads, flank_len=flank_len, np_array=False)
                if len(oneHotMat_dict['bound']) <= read_ct_thres and len(df_bound_reads):
                    oneHotMat_dict['bound'] += df_reads_2_oneHotMat(fasta, row, df_bound_reads, flank_len=flank_len, np_array=False)
                if len(oneHotMat_dict['recently bound']) <= read_ct_thres and len(df_recently_bound_reads):
                    oneHotMat_dict['recently bound'] += df_reads_2_oneHotMat(fasta, row, df_recently_bound_reads, flank_len=flank_len, np_array=False)

        if (len(oneHotMat_dict['recently bound']) >= read_ct_thres) and (len(oneHotMat_dict['bound']) >= read_ct_thres) and (len(oneHotMat_dict['unbound']) >= read_ct_thres):
            print("All labels exceed read count threshold.")
            break
    
        motif_ct += 1
        if motif_ct % 1000 == 0:
            print(f"Unbound reads: {len(oneHotMat_dict['unbound'])}, bound reads: {len(oneHotMat_dict['bound'])}, recently bound reads: {len(oneHotMat_dict['recently bound'])}")

    for edit_type in oneHotMat_dict:
        oneHotMat_dict[edit_type] = np.asarray(oneHotMat_dict[edit_type])
        print(f'{edit_type}: {oneHotMat_dict[edit_type].shape}')
    
    return(oneHotMat_dict)

def motif_row_2_df_reads_2_oneHotMat_dict_v2(bam, fasta, row, 
                                             center_trough_pos_list, left_peak_pos_list, right_peak_pos_list, 
                                             trough_thres, peak_thres, global_thres,
                                             flank_len=150, center_trough_covs_frac_thres = 0.5):

    peak_pos_list = left_peak_pos_list + right_peak_pos_list
    oneHotMat_dict = {'recently bound': [], 'bound': [], 'unbound': []}

    df_reads = motif_row_2_df_reads(bam, row)
    if len(df_reads):
        df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

        if len(df_reads):
            df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
            df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
            df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)
            df_reads = df_reads.dropna(subset=['trough_edits_frac', 'peak_edits_frac', 'global_edits_frac'])

            # label unbound bound, and recently bound reads
            df_unbound_reads = df_reads.loc[(df_reads['peak_edits_frac'] < peak_thres) & (df_reads['global_edits_frac'] < global_thres)].copy()
            df_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (df_reads['trough_edits_frac'] < trough_thres)].copy()
            df_recently_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (~df_reads.index.isin(df_bound_reads.index))].copy()

            # calculate one hot matrix for each read type
            if len(df_unbound_reads) > 0:
                oneHotMat_dict['unbound'] += df_reads_2_oneHotMat(fasta, row, df_unbound_reads, flank_len=flank_len, np_array=False)
            if len(df_bound_reads):
                oneHotMat_dict['bound'] += df_reads_2_oneHotMat(fasta, row, df_bound_reads, flank_len=flank_len, np_array=False)
            if len(df_recently_bound_reads):
                oneHotMat_dict['recently bound'] += df_reads_2_oneHotMat(fasta, row, df_recently_bound_reads, flank_len=flank_len, np_array=False)

    for edit_type in oneHotMat_dict:
        oneHotMat_dict[edit_type] = np.asarray(oneHotMat_dict[edit_type])

    return(oneHotMat_dict)


def get_read_type_frac_for_param_combos(bam, row, df_motif_features,
                                        footprint_weight_list=np.round(np.arange(-2, 2.1, 0.1), 1), peak_weight_list=np.round(np.arange(-2, 2.1, 0.1), 1), global_std_weight_list=[2],
                                        flank_len=150, center_trough_covs_frac_thres = 0.5):

    center_trough_pos_list = df_motif_features['center_pos_list'][0]
    left_peak_pos_list = df_motif_features['l_peak_pos_list'][0]
    right_peak_pos_list = df_motif_features['r_peak_pos_list'][0]
    peak_pos_list = left_peak_pos_list + right_peak_pos_list
    read_type_list = ['unbound', 'bound', 'recently bound']

    df_reads = motif_row_2_df_reads(bam, row)
    df_read_type_frac_params_combo_all = pd.DataFrame()
    if len(df_reads):
        df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

        if len(df_reads):
            df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
            df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
            df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)
            df_reads = df_reads.dropna(subset=['trough_edits_frac', 'peak_edits_frac', 'global_edits_frac'])
            
            for global_std_w in global_std_weight_list:
                global_thres = df_motif_features['inactive_mean'][0] + df_motif_features['inactive_std'][0]*global_std_w

                for peak_w in peak_weight_list:
                        peak_thres = df_motif_features['left_peak'][0]*peak_w*0.5 + df_motif_features['right_peak'][0]*peak_w*0.5 + df_motif_features['footprint'][0]*(1-peak_w)
                        
                        for footprint_w in footprint_weight_list:
                            trough_thres = df_motif_features['left_peak'][0]*footprint_w*0.5 + df_motif_features['right_peak'][0]*footprint_w*0.5 + df_motif_features['footprint'][0]*(1-footprint_w)
                            
                            # label unbound bound, and recently bound reads
                            df_unbound_reads = df_reads.loc[(df_reads['peak_edits_frac'] < peak_thres) & (df_reads['global_edits_frac'] < global_thres)].copy()
                            df_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (df_reads['trough_edits_frac'] < trough_thres)].copy()
                            df_recently_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (~df_reads.index.isin(df_bound_reads.index))].copy()

                            df_read_type_frac_params_combo = pd.DataFrame({'read_type': read_type_list, 'read_ct_defined': [len(df_unbound_reads), len(df_bound_reads), len(df_recently_bound_reads)]})
                            df_read_type_frac_params_combo['read_frac_defined'] = df_read_type_frac_params_combo['read_ct_defined'] / df_read_type_frac_params_combo['read_ct_defined'].sum()
                            df_read_type_frac_params_combo['footprint_weight'] = footprint_w
                            df_read_type_frac_params_combo['peak_weight'] = peak_w
                            df_read_type_frac_params_combo['global_std_weight'] = global_std_w
                            df_read_type_frac_params_combo_all = pd.concat([df_read_type_frac_params_combo_all, df_read_type_frac_params_combo], ignore_index=True)

    return df_read_type_frac_params_combo_all

def row_2_defined_labels_2_df_read_type_frac(bam, row, row_features,
                                             flank_len=150, center_trough_covs_frac_thres = 0.5):

    center_trough_pos_list = row_features['center_pos_list']
    left_peak_pos_list = row_features['l_peak_pos_list']
    right_peak_pos_list = row_features['r_peak_pos_list']
    peak_pos_list = left_peak_pos_list + right_peak_pos_list
    read_type_list = ['unbound', 'bound', 'recently bound']

    peak_w = row_features['peak_weight']
    footprint_w = row_features['footprint_weight']
    global_std_w = row_features['global_std_weight']
    peak_thres = row_features['left_peak']*peak_w*0.5 + row_features['right_peak']*peak_w*0.5 + row_features['footprint']*(1-peak_w)
    trough_thres = row_features['left_peak']*footprint_w*0.5 + row_features['right_peak']*footprint_w*0.5 + row_features['footprint']*(1-footprint_w)
    global_thres = row_features['inactive_mean'] + row_features['inactive_std']*global_std_w

    df_reads = motif_row_2_df_reads(bam, row)
    df_read_type_frac = pd.DataFrame()
    if len(df_reads):
        df_reads = prepare_df_reads_3readType(df_reads, row, flank_len, center_trough_pos_list, center_trough_covs_frac_thres)

        if len(df_reads):
            df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
            df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
            df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)
            df_reads = df_reads.dropna(subset=['trough_edits_frac', 'peak_edits_frac', 'global_edits_frac'])

            # label unbound bound, and recently bound reads
            df_unbound_reads = df_reads.loc[(df_reads['peak_edits_frac'] < peak_thres) & (df_reads['global_edits_frac'] < global_thres)].copy()
            df_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (df_reads['trough_edits_frac'] < trough_thres)].copy()
            df_recently_bound_reads = df_reads.loc[(~df_reads.index.isin(df_unbound_reads.index)) & (~df_reads.index.isin(df_bound_reads.index))].copy()

            df_read_type_frac = pd.DataFrame({'TF': row['TF'], 'chr': row['chr'], 'start': row['start'], 'end': row['end'], 'strand': row['strand'],
                                              'read_type': read_type_list, 'read_ct': [len(df_unbound_reads), len(df_bound_reads), len(df_recently_bound_reads)]})
            df_read_type_frac['read_frac'] = df_read_type_frac['read_ct'] / df_read_type_frac['read_ct'].sum()

    return df_read_type_frac

def prepare_df_reads_3readType_v2(df_reads, row, row_features, flank_len = 150, center_trough_covs_frac_thres = 0.5):
    df_reads = df_reads.copy()
    center_trough_pos_list = row_features['center_pos_list']
    left_peak_pos_list = row_features['l_peak_pos_list']
    right_peak_pos_list = row_features['r_peak_pos_list']
    peak_pos_list = left_peak_pos_list + right_peak_pos_list

    chrom = row['chr']
    TF_center = row['center']
    strand = row['strand']
    df_reads['chr'] = chrom
    df_reads['center'] = TF_center
    df_reads['strand'] = strand

    # get edits list
    df_reads['C2T_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['C2T']), axis=1)
    df_reads['G2A_edits'] = df_reads.apply(lambda row2: read2edits(row2['reads'], TF_center, flank_len, strand, edit_types=['G2A']), axis=1)
    df_reads['read_edit_type'] = df_reads.apply(lambda row2: get_edit_type(sum(row2['C2T_edits']), sum(row2['G2A_edits'])), axis=1)
    df_reads['edits'] = df_reads.apply(lambda row2: row2['C2T_edits'] if row2['read_edit_type']=='C2T_dominant' else row2['G2A_edits'] if row2['read_edit_type']=='G2A_dominant' else [x + y for x, y in zip(row2['C2T_edits'], row2['G2A_edits'])], axis=1)
    
    # get coverage list
    df_reads['covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand), axis=1)
    df_reads['editable_C_covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand, editable_pos_only=True, editable_bases=['C']), axis=1)
    df_reads['editable_G_covs'] = df_reads.apply(lambda row2: read2covs(row2['reads'], row2['center'], flank_len, strand, editable_pos_only=True, editable_bases=['G']), axis=1)
    df_reads['editable_covs'] = df_reads.apply(lambda row2: row2['editable_C_covs'] if row2['read_edit_type']=='C2T_dominant' else row2['editable_G_covs'] if row2['read_edit_type']=='G2A_dominant' else [x + y for x, y in zip(row2['editable_C_covs'], row2['editable_G_covs'])], axis=1)

    # filter reads by coverage in feature regions
    df_reads['center_trough_covs_frac'] = df_reads['covs'].apply(lambda x: covs_frac_in_range(x, center_trough_pos_list, flank_len))
    df_reads = df_reads.loc[(df_reads['center_trough_covs_frac']>=center_trough_covs_frac_thres)]

    # calculate feature region edit fractions
    df_reads['trough_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], center_trough_pos_list, flank_len), axis=1)
    df_reads['peak_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], peak_pos_list, flank_len), axis=1)
    df_reads['global_edits_frac'] = df_reads.apply(lambda row2: edits_frac_in_range(row2['edits'], row2['editable_covs'], [x for x in range(-flank_len, flank_len+1)], flank_len), axis=1)
    df_reads = df_reads.dropna(subset=['trough_edits_frac', 'peak_edits_frac', 'global_edits_frac']).copy()

    return(df_reads)

def df_reads_add_defined_labels(df_reads, row_features, flank_len=150):
    peak_w = row_features['peak_weight']
    footprint_w = row_features['footprint_weight']
    global_std_w = row_features['global_std_weight']
    peak_thres = row_features['left_peak']*peak_w*0.5 + row_features['right_peak']*peak_w*0.5 + row_features['footprint']*(1-peak_w)
    trough_thres = row_features['left_peak']*footprint_w*0.5 + row_features['right_peak']*footprint_w*0.5 + row_features['footprint']*(1-footprint_w)
    global_thres = row_features['inactive_mean'] + row_features['inactive_std']*global_std_w

    # label unbound bound, and recently bound reads
    df_reads = df_reads.copy()
    unbound_index = df_reads.loc[(df_reads['peak_edits_frac'] < peak_thres) & (df_reads['global_edits_frac'] < global_thres)].index
    bound_index = df_reads.loc[(~df_reads.index.isin(unbound_index)) & (df_reads['trough_edits_frac'] < trough_thres)].index
    recently_bound_index = df_reads.loc[(~df_reads.index.isin(unbound_index)) & (~df_reads.index.isin(bound_index))].index
    df_reads['read_type'] = 'unbound'
    df_reads.loc[bound_index, 'read_type'] = 'bound'
    df_reads.loc[recently_bound_index, 'read_type'] = 'recently bound'

    return df_reads


def df_loci_pair_2_df_cobinding_defined(bam, df_loci_pairs, df_motif_stats, flank_len=150, shared_read_ct_thres=100):
    read_type_list = ['unbound', 'bound', 'recently bound']
    selected_colnames = ['TF', 'chr', 'center', 'strand', 'flank_start', 'flank_end']

    # co-bdining state template
    read_type_1 = []
    read_type_2 = []
    for i in read_type_list:
        for j in read_type_list:
            read_type_1.append(i)
            read_type_2.append(j)
    df_cobinding_states_template = pd.DataFrame({'read_type_1': read_type_1, 'read_type_2': read_type_2})

    df_cobinding_all = pd.DataFrame()
    for index, row in df_loci_pairs.iterrows():
        row1 = row[[x+'1' for x in selected_colnames]]
        row2 = row[[x+'2' for x in selected_colnames]]
        row1.index = selected_colnames
        row2.index = selected_colnames

        if (row1['TF'] in df_motif_stats['TF'].tolist()) and (row2['TF'] in df_motif_stats['TF'].tolist()):
            # TF1 reads
            row_TF1 = df_motif_stats.loc[df_motif_stats['TF'] == row1['TF']].iloc[0]
            df_reads1 = motif_row_2_df_reads(bam, row1)
            if len(df_reads1):
                df_reads1 = prepare_df_reads_3readType_v2(df_reads1, row1, row_TF1, flank_len)

            # TF2 reads
            row_TF2 = df_motif_stats.loc[df_motif_stats['TF'] == row2['TF']].iloc[0]
            df_reads2 = motif_row_2_df_reads(bam, row2)
            if len(df_reads2):
                df_reads2 = prepare_df_reads_3readType_v2(df_reads2, row2, row_TF2, flank_len)

            # predict read types on shared reads
            df_reads1_shared = df_reads1.loc[df_reads1['names'].isin(df_reads2['names'])].copy()
            df_reads2_shared = df_reads2.loc[df_reads2['names'].isin(df_reads1['names'])].copy()

            if len(df_reads1_shared) < shared_read_ct_thres or len(df_reads2_shared) < shared_read_ct_thres:
                continue
            
            df_reads1_shared = df_reads_add_defined_labels(df_reads1_shared, row_TF1, flank_len)
            df_reads2_shared = df_reads_add_defined_labels(df_reads2_shared, row_TF2, flank_len)
            df_predicted_labels = pd.merge(df_reads1_shared[['names', 'read_type']], df_reads2_shared[['names', 'read_type']], on='names', suffixes=['_1', '_2'])
            row['shared_read_ct'] = len(df_reads1_shared)
            
            # expected read fractions
            uniq_value, value_ct = np.unique(df_predicted_labels['read_type_1'], return_counts=True)
            df_temp1 = pd.DataFrame({'read_type_1': uniq_value, 'read_frac_1': value_ct/np.sum(value_ct)})
            uniq_value, value_ct = np.unique(df_predicted_labels['read_type_2'], return_counts=True)
            df_temp2 = pd.DataFrame({'read_type_2': uniq_value, 'read_frac_2': value_ct/np.sum(value_ct)})
            df_cobinding_expected = pd.merge(df_cobinding_states_template, df_temp1, on='read_type_1', how='left').fillna(0)
            df_cobinding_expected = pd.merge(df_cobinding_expected, df_temp2, on='read_type_2', how='left').fillna(0)
            df_cobinding_expected['expected_read_fraction'] = df_cobinding_expected['read_frac_1'] * df_cobinding_expected['read_frac_2']

            # observed read fractions
            df_cobinding_observed = df_predicted_labels.groupby(['read_type_1', 'read_type_2']).size().reset_index().rename(columns={0:'observed_read_ct'})
            df_cobinding_observed['observed_read_fraction'] = df_cobinding_observed['observed_read_ct'] / df_cobinding_observed['observed_read_ct'].sum()

            df_cobinding = pd.merge(df_cobinding_expected, df_cobinding_observed, how='left').fillna(0)
            df_cobinding = pd.concat([pd.DataFrame([row] * len(df_cobinding)).reset_index(drop=True), df_cobinding], axis=1)
            df_cobinding_all = pd.concat([df_cobinding_all, df_cobinding])

    return df_cobinding_all

def df_loci_pair_2_df_cobinding_predicted(bam, fasta, df_loci_pairs, df_motif_stats, models_dict_all, flank_len=150, shared_read_ct_thres=100):
    read_type_list = ['unbound', 'bound', 'recently bound']
    selected_colnames = ['TF', 'chr', 'center', 'strand', 'flank_start', 'flank_end']

    # co-binding state template
    read_type_1 = []
    read_type_2 = []
    for i in read_type_list:
        for j in read_type_list:
            read_type_1.append(i)
            read_type_2.append(j)
    df_cobinding_states_template = pd.DataFrame({'read_type_1': read_type_1, 'read_type_2': read_type_2})

    df_cobinding_all = pd.DataFrame()
    for index, row in tqdm(df_loci_pairs.iterrows(), total=len(df_loci_pairs), desc='Processing loci pairs ...'):
        row1 = row[[x+'1' for x in selected_colnames]]
        row2 = row[[x+'2' for x in selected_colnames]]
        row1.index = selected_colnames
        row2.index = selected_colnames

        if (row1['TF'] in df_motif_stats['TF'].tolist()) and (row2['TF'] in df_motif_stats['TF'].tolist()):
            # TF1 reads
            row_TF1 = df_motif_stats.loc[df_motif_stats['TF'] == row1['TF']].iloc[0]
            df_reads1 = motif_row_2_df_reads(bam, row1)
            if len(df_reads1):
                df_reads1 = prepare_df_reads_3readType_v2(df_reads1, row1, row_TF1, flank_len)

            # TF2 reads
            row_TF2 = df_motif_stats.loc[df_motif_stats['TF'] == row2['TF']].iloc[0]
            df_reads2 = motif_row_2_df_reads(bam, row2)
            if len(df_reads2):
                df_reads2 = prepare_df_reads_3readType_v2(df_reads2, row2, row_TF2, flank_len)

            # predict read types on shared reads
            df_reads1_shared = df_reads1.loc[df_reads1['names'].isin(df_reads2['names'])].copy()
            df_reads2_shared = df_reads2.loc[df_reads2['names'].isin(df_reads1['names'])].copy()

            if len(df_reads1_shared) < shared_read_ct_thres or len(df_reads2_shared) < shared_read_ct_thres:
                continue

            df_probs1 = df_reads_2_predicted_probs(df_reads1_shared, model=models_dict_all[row1['TF']]['all_channel'], fasta=fasta, row=row1, flank_len=flank_len)
            df_probs2 = df_reads_2_predicted_probs(df_reads2_shared, model=models_dict_all[row2['TF']]['all_channel'], fasta=fasta, row=row2, flank_len=flank_len)
            df_predicted = pd.merge(df_probs1, df_probs2, on='names', suffixes=['_1', '_2'])
            
            row['shared_read_ct'] = len(df_reads1_shared)
            
            # expected read fractions
            uniq_value, value_ct = np.unique(df_predicted['read_type_1'], return_counts=True)
            df_temp1 = pd.DataFrame({'read_type_1': uniq_value, 'read_frac_1': value_ct/np.sum(value_ct)})
            df_temp1_prob = pd.DataFrame({'read_type_1': read_type_list, 'read_prob_1': df_probs1[['read_prob ' + read_type for read_type in read_type_list]].mean()})
            uniq_value, value_ct = np.unique(df_predicted['read_type_2'], return_counts=True)
            df_temp2 = pd.DataFrame({'read_type_2': uniq_value, 'read_frac_2': value_ct/np.sum(value_ct)})
            df_temp2_prob = pd.DataFrame({'read_type_2': read_type_list, 'read_prob_2': df_probs2[['read_prob ' + read_type for read_type in read_type_list]].mean()})
            df_cobinding_expected = pd.merge(df_cobinding_states_template, df_temp1, on='read_type_1', how='left').fillna(0)
            df_cobinding_expected = pd.merge(df_cobinding_expected, df_temp1_prob, on='read_type_1', how='left')
            df_cobinding_expected = pd.merge(df_cobinding_expected, df_temp2, on='read_type_2', how='left').fillna(0)
            df_cobinding_expected = pd.merge(df_cobinding_expected, df_temp2_prob, on='read_type_2', how='left')
            
            df_cobinding_expected['expected_read_fraction'] = df_cobinding_expected['read_frac_1'] * df_cobinding_expected['read_frac_2']
            df_cobinding_expected['expected_probability'] = df_cobinding_expected['read_prob_1'] * df_cobinding_expected['read_prob_2']
            
            # observed read fractions
            df_cobinding_observed = df_predicted.groupby(['read_type_1', 'read_type_2']).size().reset_index().rename(columns={0:'observed_read_ct'})
            df_cobinding_observed['observed_read_fraction'] = df_cobinding_observed['observed_read_ct'] / df_cobinding_observed['observed_read_ct'].sum()
            df_cobinding_observed = pd.merge(df_cobinding_states_template, df_cobinding_observed, on=['read_type_1', 'read_type_2'], how='left').fillna(0)
            
            observed_probs = []
            for index, row_obs in df_cobinding_observed.iterrows():
                read_type_1 = row_obs['read_type_1']
                read_type_2 = row_obs['read_type_2']
                temp = (df_predicted['read_prob ' + read_type_1 + '_1'] * df_predicted['read_prob ' + read_type_2 + '_2']).mean()
                observed_probs.append(temp)
            df_cobinding_observed['observed_probability'] = observed_probs

            df_cobinding = pd.merge(df_cobinding_expected, df_cobinding_observed, on=['read_type_1', 'read_type_2'], how='left').fillna(0)
            df_cobinding = pd.concat([pd.DataFrame([row] * len(df_cobinding)).reset_index(drop=True), df_cobinding], axis=1)
            df_cobinding_all = pd.concat([df_cobinding_all, df_cobinding])

    return df_cobinding_all

def adjust_motif_center_v2(row, row_features):
    # center_adj_TF = (max(row_features['center_pos_list']) + min(row_features['center_pos_list']))/2
    center_adj_TF = np.mean(row_features['center_pos_list'])
    # center_adj_TF = np.median(row_features['center_pos_list'])
    center_TF = row['center']
    center_TF = (center_TF + center_adj_TF) if (row['strand'] == '+') else (center_TF- center_adj_TF)
    return center_TF

def calc_center_distance(row_loci_pairs, row_TF1, row_TF2):
    center1 = row_loci_pairs['center1']
    strand1 = row_loci_pairs['strand1']
    # center1_adj = (max(row_TF1['center_pos_list']) + min(row_TF1['center_pos_list']))/2
    center1_adj = np.mean(row_TF1['center_pos_list'])
    # center1_adj = np.median(row_TF1['center_pos_list'])
    center1_span = max(row_TF1['center_pos_list']) - min(row_TF1['center_pos_list'])
    if strand1 == '+':
        center1_corrected = center1 + center1_adj
    else:
        center1_corrected = center1 - center1_adj

    center2 = row_loci_pairs['center2']
    strand2 = row_loci_pairs['strand2']
    # center2_adj = (max(row_TF2['center_pos_list']) + min(row_TF2['center_pos_list']))/2
    center2_adj = np.mean(row_TF2['center_pos_list'])
    # center2_adj = np.median(row_TF2['center_pos_list'])
    center2_span = max(row_TF2['center_pos_list']) - min(row_TF2['center_pos_list'])
    if strand2 == '+':
        center2_corrected = center2 + center2_adj
    else:
        center2_corrected = center2 - center2_adj

    if abs(center1_corrected - center2_corrected) > (center1_span + center2_span)/2:
        if strand1 == '+':
            return center2_corrected - center1_corrected
        else:
            return center1_corrected - center2_corrected
    else:
        return np.nan


def calc_center_distance_v2(row_loci_pairs, df_motif_stats, adj_method_num):
    row_TF1 = df_motif_stats.loc[df_motif_stats['TF'] == row_loci_pairs['TF1']].iloc[0]
    center1 = row_loci_pairs['center1']
    strand1 = row_loci_pairs['strand1']
    if adj_method_num == 0:
        center1_adj = (max(row_TF1['center_pos_list']) + min(row_TF1['center_pos_list']))/2
    if adj_method_num == 1:
        center1_adj = (max(row_TF1['center_pos_list']) + min(row_TF1['center_pos_list']) + max(row_TF1['l_peak_pos_list']) + min(row_TF1['l_peak_pos_list']) + max(row_TF1['r_peak_pos_list']) + min(row_TF1['r_peak_pos_list']))/6
    if adj_method_num == 2:
        center1_adj = np.median(row_TF1['center_pos_list'])
    if adj_method_num == 3:
        center1_adj = np.mean(row_TF1['center_pos_list'])
    if adj_method_num == 4:
        center1_adj = np.median(row_TF1['center_pos_list'] + row_TF1['l_peak_pos_list'] + row_TF1['r_peak_pos_list'])
    if adj_method_num == 5:
        center1_adj = np.mean(row_TF1['center_pos_list'] + row_TF1['l_peak_pos_list'] + row_TF1['r_peak_pos_list'])

    center1_span = max(row_TF1['center_pos_list']) - min(row_TF1['center_pos_list'])
    if strand1 == '+':
        center1_corrected = center1 + center1_adj
    else:
        center1_corrected = center1 - center1_adj

    row_TF2 = df_motif_stats.loc[df_motif_stats['TF'] == row_loci_pairs['TF2']].iloc[0]
    center2 = row_loci_pairs['center2']
    strand2 = row_loci_pairs['strand2']
    if adj_method_num == 0:
        center2_adj = (max(row_TF2['center_pos_list']) + min(row_TF2['center_pos_list']))/2
    if adj_method_num == 1:
        center2_adj = (max(row_TF1['center_pos_list']) + min(row_TF1['center_pos_list']) + max(row_TF1['l_peak_pos_list']) + min(row_TF1['l_peak_pos_list']) + max(row_TF1['r_peak_pos_list']) + min(row_TF1['r_peak_pos_list']))/6
    if adj_method_num == 2:
        center2_adj = np.median(row_TF2['center_pos_list'])
    if adj_method_num == 3:
        center2_adj = np.mean(row_TF2['center_pos_list'])
    if adj_method_num == 4:
        center2_adj = np.median(row_TF1['center_pos_list'] + row_TF1['l_peak_pos_list'] + row_TF1['r_peak_pos_list'])
    if adj_method_num == 5:
        center2_adj = np.mean(row_TF1['center_pos_list'] + row_TF1['l_peak_pos_list'] + row_TF1['r_peak_pos_list'])

    center2_span = max(row_TF2['center_pos_list']) - min(row_TF2['center_pos_list'])
    if strand2 == '+':
        center2_corrected = center2 + center2_adj
    else:
        center2_corrected = center2 - center2_adj

    if abs(center1_corrected - center2_corrected) > (center1_span + center2_span)/2:
        if strand1 == '+':
            return center2_corrected - center1_corrected
        else:
            return center1_corrected - center2_corrected
    else:
        return np.nan
    
def calc_fft(x, y, randomization_test=False, num_iterations=1000, alpha=0.05):
    # Remove NaNs and infinite values
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite_mask], y[finite_mask]

    # Perform FFT
    yf = fft(y)
    xf = np.linspace(0.0, 1.0 / (2.0 * x.diff().mean()), len(y) // 2)
    fft_magnitude = 2.0 / len(y) * np.abs(yf[:len(y) // 2])

    # Create DataFrame for FFT results
    df_fft = pd.DataFrame({
        'frequency': xf,
        'amplitude': fft_magnitude
    })
    df_fft = df_fft.loc[df_fft['frequency'] > 0].copy()
    df_fft['period'] = 1 / df_fft['frequency']

    # Perform Randomization Test if specified
    if randomization_test:
        # Initialize array to store shuffled amplitudes
        randomized_amplitudes = np.zeros((num_iterations, len(df_fft)))
        
        # Randomization (Shuffling) loop
        for i in range(num_iterations):
            y_shuffled = np.random.permutation(y)
            yf_shuffled = fft(y_shuffled)
            fft_magnitude_shuffled = 2.0 / len(y_shuffled) * np.abs(yf_shuffled[:len(y) // 2])
            randomized_amplitudes[i, :] = fft_magnitude_shuffled[:len(df_fft)]

        # Calculate p-values and the 95th percentile threshold
        df_fft['pval'] = np.mean(randomized_amplitudes >= df_fft['amplitude'].values, axis=0)
        df_fft['CI_thres'] = np.percentile(randomized_amplitudes, 100 * (1 - alpha), axis=0)

    return df_fft

def calc_fft(x, y, randomization_test=False, num_iterations=1000, alpha=0.05):
    # Remove NaNs and infinite values
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite_mask], y[finite_mask]

    # Perform FFT
    yf = fft(y)
    xf = np.linspace(0.0, 1.0 / (2.0 * x.diff().mean()), len(y) // 2)[1:]  # Exclude frequency 0
    fft_magnitude = (2.0 / len(y) * np.abs(yf[:len(y) // 2]))[1:]  # Exclude amplitude at frequency 0

    # Create DataFrame for FFT results, excluding frequency 0
    df_fft = pd.DataFrame({
        'frequency': xf,
        'amplitude': fft_magnitude
    })
    df_fft['period'] = 1 / df_fft['frequency']

    # Perform Randomization Test if specified
    if randomization_test:
        # Initialize array to store shuffled amplitudes
        randomized_amplitudes = np.zeros((num_iterations, len(df_fft)))
        
        # Randomization (Shuffling) loop
        for i in range(num_iterations):
            y_shuffled = np.random.permutation(y)
            yf_shuffled = fft(y_shuffled)
            fft_magnitude_shuffled = (2.0 / len(y_shuffled) * np.abs(yf_shuffled[:len(y) // 2]))[1:]  # Exclude freq 0
            randomized_amplitudes[i, :] = fft_magnitude_shuffled[:len(df_fft)]

        # Calculate p-values and the 95th percentile threshold
        df_fft['pval'] = np.mean(randomized_amplitudes >= df_fft['amplitude'].values, axis=0)
        df_fft['CI_thres'] = np.percentile(randomized_amplitudes, 100 * (1 - alpha), axis=0)

    return df_fft


def df_cobinding_2_df_fft(df_cobinding, x_column='center_dist', y_column='delta_oe_by_pos', randomization_test=False, num_iterations=1000, alpha=0.05):
    df_cobinding = df_cobinding.dropna(subset=[y_column]).copy()
    df_combined_fft = pd.DataFrame()
    for side in ['Positive', 'Negative']:
        if side == 'Positive':
            df_side = df_cobinding.loc[df_cobinding['center_dist']>0].copy()
        else:
            df_side = df_cobinding.loc[df_cobinding['center_dist']<0].copy()

        if len(df_side) == 0:
            continue
        
        df_side_fft = calc_fft(x=df_side[x_column], y=df_side[y_column].values, randomization_test=randomization_test, num_iterations=num_iterations, alpha=alpha)
        if randomization_test:
            df_side_fft['signal_noise_ratio'] = df_side_fft['amplitude'] / df_side_fft['CI_thres']
        df_side_fft['side'] = side
        df_combined_fft = pd.concat([df_combined_fft, df_side_fft], ignore_index=True)

    return df_combined_fft


def df_cobinding_2_df_fft_v2(df_cobinding, x_column='center_dist', y_column='delta_oe_prob_norm', y_column_rename='delta_oe_by_pos_norm',
                             proximity_range=[16, 100], one_side_window_size=0, 
                             randomization_test=False, num_iterations=1000, alpha=0.05, n_jobs=-1, parallel=True):
    
    max_center_dist = proximity_range[1]
    df_combined_fft = pd.DataFrame()
    for side in ['Positive', 'Negative']:
        if side == 'Positive':
            df_side = df_cobinding.loc[df_cobinding['center_dist'] > 0].copy()
        else:
            df_side = df_cobinding.loc[df_cobinding['center_dist'] < 0].copy()

        x = []
        y = []
        for center_dist in range(-max_center_dist, max_center_dist + 1):
            x.append(center_dist)
            observed_expected_delta_ma = df_side.loc[
                df_side['center_dist'].between(center_dist - one_side_window_size - 0.5, 
                                               center_dist + one_side_window_size + 0.5, inclusive='left'), 
                y_column].median()
            y.append(observed_expected_delta_ma)

        # Remove NaNs and infinite values
        x = np.array(x)
        y = np.array(y)
        finite_mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[finite_mask], y[finite_mask]

        # Perform FFT
        yf = fft(y)
        xf = np.linspace(0.0, 1.0 / (2.0 * np.diff(x).mean()), len(y) // 2)[1:]
        fft_magnitude = (2.0 / len(y) * np.abs(yf[:len(y) // 2]))[1:]  # Exclude amplitude at frequency 0

        # Create DataFrame for FFT results, excluding frequency 0
        df_side_fft = pd.DataFrame({'frequency': xf, 'amplitude': fft_magnitude})
        df_side_fft['period'] = 1 / df_side_fft['frequency']
        df_side_fft['side'] = side

        if randomization_test:
            # Define randomization logic
            def randomize_fft(i):
                df_side_shuffled = df_side.copy()
                df_side_shuffled[y_column] = np.random.permutation(df_side[y_column].values)

                x_shuffled = []
                y_shuffled = []
                for center_dist in range(-max_center_dist, max_center_dist + 1):
                    x_shuffled.append(center_dist)
                    observed_expected_delta_ma = df_side_shuffled.loc[
                        df_side_shuffled['center_dist'].between(center_dist - one_side_window_size - 0.5, 
                                                                center_dist + one_side_window_size + 0.5, inclusive='left'), 
                        y_column].median()
                    y_shuffled.append(observed_expected_delta_ma)

                # Remove NaNs and infinite values
                x_shuffled = np.array(x_shuffled)
                y_shuffled = np.array(y_shuffled)
                finite_mask = np.isfinite(x_shuffled) & np.isfinite(y_shuffled)
                x_shuffled, y_shuffled = x_shuffled[finite_mask], y_shuffled[finite_mask]

                yf_shuffled = fft(y_shuffled)
                fft_magnitude_shuffled = (2.0 / len(y_shuffled) * np.abs(yf_shuffled[:len(y_shuffled) // 2]))[1:]  # Exclude freq 0
                return fft_magnitude_shuffled[:len(df_side_fft)]

            if parallel:
                randomized_amplitudes = Parallel(n_jobs=n_jobs)(
                    delayed(randomize_fft)(i) for i in range(num_iterations)
                )
            else:
                randomized_amplitudes = np.array([
                    randomize_fft(i) for i in range(num_iterations)
                ])
            randomized_amplitudes = np.array(randomized_amplitudes)

            # Calculate p-values and the 95th percentile threshold
            df_side_fft['pval'] = np.mean(randomized_amplitudes >= df_side_fft['amplitude'].values, axis=0)
            df_side_fft['CI_thres'] = np.percentile(randomized_amplitudes, 100 * (1 - alpha), axis=0)
            df_side_fft['signal_noise_ratio'] = df_side_fft['amplitude'] / df_side_fft['CI_thres']

        df_combined_fft = pd.concat([df_combined_fft, df_side_fft], ignore_index=True)
    return df_combined_fft


def get_dominant_period(df_fft, min_period=5, max_period=50, alpha=0.05):
    side_list = ['Positive', 'Negative']
    dom_freq = []
    dom_amp = []
    signal_noise_ratio = []
    pval = []
    for side in side_list:
        df_fft_by_side = df_fft.loc[(df_fft['frequency'].between(1/max_period, 1/min_period, inclusive='both')) & 
                                    (df_fft['pval'] < alpha) &
                                    (df_fft['side'] == side)].copy()
        if len(df_fft_by_side) == 0:
            dom_freq.append(np.nan)
            dom_amp.append(np.nan)
            signal_noise_ratio.append(np.nan)
            pval.append(np.nan)
        else:
            dom_amp.append(df_fft_by_side['amplitude'].max())
            dom_freq.append(df_fft_by_side['frequency'].loc[df_fft_by_side['amplitude'].idxmax()])
            signal_noise_ratio.append(df_fft_by_side['signal_noise_ratio'].loc[df_fft_by_side['amplitude'].idxmax()])
            pval.append(df_fft_by_side['pval'].loc[df_fft_by_side['amplitude'].idxmax()])
    df_fft_stats = pd.DataFrame({'side': side_list, 'dom_freq': dom_freq, 'dom_amp': dom_amp, 'signal_noise_ratio': signal_noise_ratio, 'pval': pval})
    df_fft_stats['period'] = 1/df_fft_stats['dom_freq']
    return df_fft_stats

def get_period_shift(df_fft_stats, df_cobinding_bound_by_pos, y_column='delta_oe_by_pos'):
    min_period_shift_list = []
    df_trough_pos = pd.DataFrame()
    for index, row in df_fft_stats.iterrows():
        side = row['side']
        period = row['period']

        if 'motif_pair_orientation' in df_fft_stats.columns:
            motif_pair_orientation = row['motif_pair_orientation']
            df_cobinding_bound_by_pos2 = df_cobinding_bound_by_pos.loc[(df_cobinding_bound_by_pos['motif_pair_orientation'] == motif_pair_orientation)].copy()
        else:
            df_cobinding_bound_by_pos2 = df_cobinding_bound_by_pos.copy()

        if pd.isna(period):
            min_period_shift_list.append(np.nan)
            continue
        else:
            period_step = period / int(period+1)

        if side == 'Negative':
            df_cobinding_bound_by_pos_side = df_cobinding_bound_by_pos2.loc[(df_cobinding_bound_by_pos2['center_dist'] < 0)].copy()
        else:
            df_cobinding_bound_by_pos_side = df_cobinding_bound_by_pos2.loc[(df_cobinding_bound_by_pos2['center_dist'] > 0)].copy()
        pos_min = df_cobinding_bound_by_pos_side['center_dist'].min()
        pos_max = df_cobinding_bound_by_pos_side['center_dist'].max()

        period_shift_list = np.arange(0, period, period_step)
        mean_delta_oe_by_pos_list = []
        for period_shift in period_shift_list:
            period_shift_pos_list = [round(x) for x in np.arange(pos_min+period_shift, pos_max, period)]
            delta_values = df_cobinding_bound_by_pos_side.loc[df_cobinding_bound_by_pos_side['center_dist'].isin(period_shift_pos_list), y_column].replace([np.inf, -np.inf], np.nan).dropna()
            mean_delta_oe_by_pos_list.append(delta_values.mean())

        min_period_shift = period_shift_list[np.argmin(mean_delta_oe_by_pos_list)]
        min_period_shift_list.append(min_period_shift)

        df_trough_pos_side = pd.DataFrame({'center_dist': [round(x) for x in np.arange(pos_min+min_period_shift, pos_max, period)]})
        if 'motif_pair_orientation' in df_fft_stats.columns:
            df_trough_pos_side['motif_pair_orientation'] = motif_pair_orientation
        df_trough_pos_side['side'] = side
        df_trough_pos = pd.concat([df_trough_pos, df_trough_pos_side], ignore_index=True)
    
    df_fft_stats['period_shift'] = min_period_shift_list
    return df_fft_stats, df_trough_pos

        
