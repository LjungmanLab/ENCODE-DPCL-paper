# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 14:53:22 2021

Description
-----------
This script defines the following functions which are used to generate decreasing
fused segments file (.fused_dec.bed) from parent segments.actual_RPKM.bed file using pandas:
    
    load_from_bed: loads in bed file as pandas compatible dataframe
    split_dataframe: splits df on strand (pos or neg) for ease of merging seg on rows
    simple_merge: defines merge/fusing conditions for a simple merge of all consecutive, decreasing segments
    flexible_merge: defines merge/fusing conditions that allows one improper segment (i.e. not decreasing)
    merge_on_score: merges/fuses segments according to decreasing score and constructs new df rows
    delete_if_exists: deletes file if already exists in working directory
    cat_pos_neg: concatenates pos/neg merged/fused dfs to create final output file
    
Version Notes
-------------
v2: added requirement that chr match between segments being fused
v3: made it so segments end as soon as score hits zero (no zero scores fused into prev_seg) and have ID output correct name
v4: make cutoff of segment end score actually 1 (segs fused have score 2-9)
v5: fix bugs related to merging
v6: change cutoff of segment end score to be greater than 2 (segs fused have score 3-9); also remove commented out code
v7: *MAJOR CHANGE* Added flexible_merge and simple_merge functions, and changed method of calling merge_func
v8: *FINAL* Removed flexible_merge function (only simple will be used) and other unused funcs; fixed segment length calculation;
    changed output format to include only relevant columns and file naming convention

@author: abmcs
"""

# load data as pandas dataframe

import argparse
import os
from tqdm import tqdm
import pandas as pd
import sys


def load_from_bed(filepath):
    header = ['chr', 'start', 'end', 'ID', 'score', 'strand', 'overlapping_gene_count', 'feat_type', 'feat_pos', 'seg_length', 'pseudocounts', 'read_density', 'rpkm']
    df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    
    return df

    
def split_dataframe(df):
	pos_df = df[df.strand == "+"]
	neg_df = df[df.strand == "-"]
	return pos_df, neg_df

     
def simple_merge(row, df, prev_score, prev_chr):
    score = row.score
    _chr = row.chr
    is_dec = score <= prev_score and score > 2 and prev_chr == _chr
    return is_dec
     
     
def merge_on_score(df):
    
    '''
    Parameters
    ----------
    df : pandas dataframe
    merge_func : function that defines merging conditions
    
    Returns
    -------
    pandas df: with all merged rows
    
    Description
    -----------
    This function will take a segments file formatted as a pandas df, merge sequential rows of
    that df according to the score assigned to that segment (defined by merge_func), and output
    the following cols (where * denotes new or recalculated cols):
        chr -------------------------------- chromosome
        start ------------------------------ new segment start coordinate
        end* ------------------------------- new segment end coordinate
        ID* -------------------------------- new segment ID for fused segment
        score* ----------------------------- new segment starting score
        strand ----------------------------- strand
        feat_type* ------------------------- segment feat type (fused_dec_segment)
        seg_length* ------------------------ length of new segment (bp)
        fused_seg_ids* --------------------- list of all segments fused to create new seg
        min_segment_score* ----------------- new segment end score
        num_fused_seg* --------------------- number of segments fused to create new seg
    '''
    
    merged_rows = []
    prev_score = None
    prev_chr = None
    merge_id_col = "fused_seg_ids"
    num_seg_col = "num_fused_seg"
    end_score_col = "min_segment_score"
    seg_len_col = "seg_length"
    feat_type_col = "feat_type"
    
    if df.strand.iloc[0] == '-':
        end_col = "start"
    else:
        end_col = "end"
    
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df)):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        score = row.score
        _chr = row.chr
        seg_len = row.seg_length
        if (prev_score is not None) and simple_merge(row, df, prev_score, prev_chr):  
            prev_row = merged_rows[-1]
            prev_row[merge_id_col] += ',' + row.ID
            prev_row[end_col] = row[end_col]
            prev_row[end_score_col] = score
            prev_row[num_seg_col] += 1
            prev_row[seg_len_col] += seg_len
        else:
            row[merge_id_col] = row.ID
            row[end_score_col] = None
            merged_rows.append(row)
            row[num_seg_col] = 1
        prev_score = score
        prev_chr = _chr
    
    
    pd_df = pd.DataFrame(merged_rows)
    pd_df[feat_type_col] = "fused_dec_segment"
    pd_df[seg_len_col] = pd_df["end"] - pd_df["start"]
    
    formatted_df = pd_df[['chr', 'start', 'end', 'ID', 'score', 'strand', 'feat_type', 'seg_length', 'fused_seg_ids', 'min_segment_score', 'num_fused_seg']].rename(columns={'chr': '#chr'})
    
    return formatted_df


def delete_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    return False


def cat_pos_neg(pos_merged_df, neg_merged_df):
    seg_dfs = [pos_merged_df, neg_merged_df]
    fused_dec_df = pd.concat(seg_dfs)
    return fused_dec_df


def main(bed_filepath):
    
    '''
    Parameters
    ----------
    bed_filepath: path to segments.actual_RPKM file to be processed
    
    Description
    -----------
    Opens file, splits df on strand, reverses neg strand df order (to allow consecutive row merge), fuses segments,
    merges pos and neg dfs, and writes new fused_dec bed file.
    
    Variables
    ---------
    No variable to change this main function.
    '''
    
    
    df = load_from_bed(bed_filepath)
    	
    pos_df, neg_df = split_dataframe(df)
    
    neg_df = neg_df[::-1]
    
    pos_merged_df = merge_on_score(pos_df)
    neg_merged_df = merge_on_score(neg_df)
    
    neg_merged_df = neg_merged_df.sort_index()
    
    fused_dec_df = cat_pos_neg(pos_merged_df, neg_merged_df)
    
    fused_dec_df['ID'] = fused_dec_df['#chr']+":"+fused_dec_df['start'].astype(str)+"-"+fused_dec_df['end'].astype(str)
        
    bed_output = os.path.splitext(bed_filepath)[0] + ".fused_dec.bed"
    delete_if_exists(bed_output)
    fused_dec_df.to_csv(bed_output, sep='\t', header=True, index=False)


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fn',
                        help='Path to segments file (expects .segments.actual.bed file).')
    
    args = parser.parse_args()
        
    main(args.fn)
	
