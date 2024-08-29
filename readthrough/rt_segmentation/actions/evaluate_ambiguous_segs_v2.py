# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 09:06:55 2022

Description
-----------
This script defines the following functions which are used to... :
    
    load_files: loads in all files needed for this script
    delete_if_exists: deletes file if already exists in working directory
    evaluate_segments: evaluates and classifies previously ambiguous segments
    update_segment_coordinates: changes start/end coordinates of any previously ambiguous segments
    redefine_segment_coordinates: primary workflow for reassignment of RT_seg coordinates after reevaluation
    refine_classifications: removes redundancies, refines classifications, and formats output df

@author: abmcs
"""


import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
from bisect import bisect
import sys
import warnings


def load_files(filepath):
    
    '''
    loads in files and adds header if needed
    '''
    
    if 'eval' in filepath: # final_df containing segments that do not require eval
        # df = pd.read_csv(filepath, sep='\t', header=0, engine='python')
        df0 = pd.read_csv(filepath, sep='\t', engine='python')
        df = df0.rename(columns={'#RT_chr': 'RT_chr'})
        
    else: # eval region files
        header = ['chr', 'start', 'end', 'RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL', 'feat_length', 'counts', 'density', 'rpkm']
        
        df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
        
    return df


def delete_if_exists(filepath):
    
    '''
    removes existing output file if filename matches one in directory
    '''
    
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    
    return False


def evaluate_segments(start_df, end_df, gene_df):
    
    '''
    Parameters
    ----------
    start_df: file containing counts for the 1kb region at the start of the RT_segs requiring evaluation
    end_df: file containing counts for the 1kb region of the RT_seg before the start of the downstream gene for 
            RT_segs requiring evaluation
    gene_df: file containing counts for 1kb into the start of the downstream gene overlapped by RT_segs requiring
            evaluation
    
    Description
    -----------
    This function takes in RT_start, RT_end, and DS_gene information, calculates RT_end/RT_start and RT_end/DS_gene
    ratios, evaluates these read-through regions, and classifies each segment as either IId, IIf, or IIIa; output
    is a new df with additional cols for segments that can be reevaluated.
    
        RULES:
            1.) if E/G ratio is < 0.5, can eval, set end of RT seg as start of downstream gene
            2.) if counts at RT_end or start of gene = 0, can eval, set end of RT seg as start of downstream gene
            3.) if E/G ratio is > 0.5 and expression is not low, cannot eval, change RT_class to reflect this (IIIa)
            4.) if rpkm of gene and RT_end are really low (but not zero and don't meet E/G cutoff), just leave
                segment as is, classify as such (IIf))
        
    Variables
    ---------
    Any of the cutoffs for the evaluation rules listed above can be changed. Additionally, there is currently no
    requirement that the start of a RT_segment must be less than or equal to the end, this could be changed if 
    deemed necessary.
    '''
     
    # merge the three files
    merged_df = pd.merge(start_df, end_df, on=['RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL', 'feat_length'], how='outer', suffixes=('_start', '_end'))
    merged_df = pd.merge(merged_df, gene_df, on=['RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL', 'feat_length'], how='outer')
    merged_df.rename(columns={'chr':'chr_gene', 'start':'start_gene', 'end':'end_gene', 'counts':'counts_gene', 'density':'density_gene', 'rpkm':'rpkm_gene'}, inplace=True)
    
    # reorder columns
    merged_df = merged_df[['chr_start', 'start_start', 'end_start', 'chr_end', 'start_end', 'end_end', 'chr_gene', 'start_gene', 'end_gene',
                           'RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL', 'feat_length',
                           'counts_start', 'density_start', 'rpkm_start', 'counts_end', 'density_end', 'rpkm_end', 'counts_gene', 'density_gene', 'rpkm_gene']]
    
    # caluclate ratios of RT_end/RT_start and DS_TSS/RT_end
    merged_df['counts_E/S'] = merged_df['counts_end'] / merged_df['counts_start']
    merged_df['counts_E/G'] = merged_df['counts_end'] / merged_df['counts_gene']
    merged_df['rpkm_E/S'] = merged_df['rpkm_end'] / merged_df['rpkm_start']
    merged_df['rpkm_E/G'] = merged_df['rpkm_end'] / merged_df['rpkm_gene']
    
    # filter segments based on these ratios (need to explore these)
        
    rule_1 = merged_df[merged_df['rpkm_E/G'] < 0.5].copy()  
    rule_2 = merged_df[(merged_df['counts_end'] == 0) | (merged_df['counts_gene'] == 0)].copy()
    rule_3 = merged_df[(merged_df['rpkm_E/G'] > 0.5) & (merged_df['rpkm_end'] > 0.25) & (merged_df['rpkm_gene'] > 0.25)].copy()
    rule_4 = merged_df[(merged_df['rpkm_E/G'] > 0.5) & (merged_df['rpkm_end'] < 0.25) & (merged_df['rpkm_end'] > 0) & (merged_df['rpkm_gene'] < 0.25) & (merged_df['rpkm_gene'] > 0)].copy()
    
    # classify segments as ambiguous or not
    # merge rule 1 and 2 and create new RT_coordinates and add RT_class
    can_eval = pd.merge(rule_1, rule_2, how='outer')
    rule12_rows = []
    # for row_ind, row_vals in tqdm(enumerate(can_eval.iterrows()), total=len(can_eval)):
    for row_ind, row in can_eval.iterrows():
        # row = row_vals[1]
        row['new_chr'] = row['chr_start']
        if row['strand'] == '+':
            row['new_start'] = row['start_start']
            row['new_end'] = row['start_gene']
            row['new_ID'] = row['new_chr']+':'+str(row['new_start'])+'-'+str(row['new_end'])
            rule12_rows.append(row)
        else:
            row['new_start'] = row['end_gene']
            row['new_end'] = row['end_start']
            row['new_ID'] = row['new_chr']+':'+str(row['new_start'])+'-'+str(row['new_end'])
            rule12_rows.append(row)
    
    rule12_df = pd.DataFrame(rule12_rows)
    rule12_df['RT_class'] = 'IId'
    rule12_df['new_length'] = rule12_df['new_end'] - rule12_df['new_start']
    
    # add RT_class to rule 3
    rule_3['RT_class'] = 'IIIa'
    
    # add RT_class to rule 4
    rule_4['RT_class'] = 'IIf'
    
    # cat rule 3 and 4 and add "new" coordinate columns (same as previous RT_seg)
    ### can be some overlap of RT_IDs
    rule34_df = pd.concat([rule_3, rule_4], ignore_index=True)
    rule34_df[['new_chr', 'intermediate']] = rule34_df['RT_ID'].str.split(':', n=1, expand=True)
    rule34_df[['new_start', 'new_end']] = rule34_df['intermediate'].str.split('-', n=1, expand=True)
    rule34_df.drop('intermediate', inplace=True, axis=1)
    rule34_df['new_ID'] = rule34_df['RT_ID']
    rule34_df['new_length'] = np.nan
    
    # cat all rules
    fixed_df = pd.concat([rule12_df, rule34_df], ignore_index=True)
    fixed_df = fixed_df.astype({'new_start': int, 'new_end': int})
    fixed_df = fixed_df.sort_values(by = ['new_chr', 'new_start'])
    output_df = fixed_df[['new_chr', 'new_start', 'new_end', 'new_ID', 'score', 'strand', 'new_length', 'RT_ID', 'RT_class', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL']]
    
    return output_df


def update_segment_coordinates(mult_df, index_list):
    '''
    Parameters
    ----------
    mult_df: a dataframe containing all rows with matching RT_IDs (must be greater than 1 row), rows have been
            ordered by coordinates (so position in df matches order of genes)
    index_list: list of indeces in mult_df that have been assigned class IId (reevaluated)

    Description
    -------
    This function takes a dataframe of all rows for each unique RT_ID with multiple gene overlaps that contains a
    segment of class IId (can be reevaluated) and redefines the segment coordinates of that segment and redefines
    the coordinates of any additional segments that are created for overlapping genes from the original segment.
    Comments describe logic being used for coordinate reassignment.
    
    Variables
    -------
    There are instances where previous cutoffs are used in this function (gene length cutoff = 1000, overlapping
    gene cutoff = 1000, and so on). If any changes are made upstream, they will need to be modified here as well.
    '''
    
    for idx, row in mult_df.iterrows():
        # check if index is in index_list
        if idx in index_list:
            continue
        # if idx is not in index_list, these are segments that need to be changed based on the IId segments
        else:
            position = bisect(index_list, idx)
            # if bisect = 0, then it's the first ol_gene
            # start stays the same but end changes to start of first ds gene that is was evaled
            # there should be no added rows 'X' before a 'IId' reeval segment
            if position == 0:
                before_idx = None
                after_idx = index_list[position]
                # if small or not_expressed, set to ds reeval segment
                if (mult_df.loc[idx, 'expressed'] == 'no') | (mult_df.loc[idx, 'ol_gene_length'] < 1000):
                    mult_df.loc[idx, 'new_chr'] = mult_df.loc[after_idx, 'RT_chr']
                    mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'new_start']
                    mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'new_end']
                    mult_df.loc[idx, 'new_ID'] = mult_df.loc[after_idx, 'new_ID']
                    mult_df.loc[idx, 'score'] = mult_df.loc[after_idx, 'score']
                    mult_df.loc[idx, 'strand'] = mult_df.loc[after_idx, 'strand']
                    mult_df.loc[idx, 'new_length'] = mult_df.loc[after_idx, 'new_length']
                    # make sure RT_length is a positive value (would not be if segment ends in middle of ds gene)
                    if mult_df.loc[idx, 'new_length'] < 0:
                        mult_df = mult_df.drop([idx])
                        continue
                    # overwrite original RT_classes if small or not expressed (these will be checked and reassigned later)
                    if (mult_df.loc[idx, 'ol_gene_length'] < 1000):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIa' + ',' + mult_df.loc[after_idx, 'RT_class_y']
                        continue
                    if (mult_df.loc[idx, 'expressed'] == 'no'):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIc' + ',' + mult_df.loc[after_idx, 'RT_class_y']
                        continue
                # if no entries and not small/not expressed
                if pd.isna(mult_df.loc[idx, 'new_start']):
                    # change stranded column values
                    if mult_df.loc[idx, 'RT_strand'] == '+':
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[idx, 'RT_start']
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                    else:
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[idx, 'RT_end']
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                    # change non-stranded column values
                    mult_df.loc[idx, 'new_chr'] = mult_df.loc[idx, 'RT_chr']
                    mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                    mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    mult_df.loc[idx, 'num_geneOL'] = len(mult_df) - after_idx
                    # make sure RT_length is a positive value (would not be if segment ends in middle of ds gene)
                    if mult_df.loc[idx, 'new_length'] < 0:
                        mult_df = mult_df.drop([idx])
                        continue
                    # add RT_class as IId, because this was changed as a result of reeval
                    mult_df.loc[idx, 'RT_class_y'] = 'IId'
                # row could be IIf or IIIa; end can change, but start will stay the same
                else:
                    # if IIf is in this position, it is incorrect (met criteria but not partial overlap), change
                    # to IIIa and change coordinates
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIf':
                        warnings.warn("Exception: segment class IIf in first position. Class changed to IIIa.")
                        mult_df.loc[idx, 'RT_class_y'] = 'IIIa'
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                        else:
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    # if class is IIIa, change coordinates
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIIa':
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                        else:
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
            # if bisect = len(index_list), it's the last ol_gene;
            # start changes to the end of the nearest ds gene that was evaled and end stays the same
            if position == len(index_list):
                before_idx = index_list[position - 1]
                after_idx = None
                if mult_df.loc[idx, 'RT_class_y'] == 'X':
                    for i, r in mult_df.iloc[idx+1:len(mult_df)].iterrows():
                        if pd.isna(mult_df.loc[i, 'RT_class_y']):
                            mult_df.loc[i, 'new_chr'] = mult_df.loc[idx, 'new_chr']
                            mult_df.loc[i, 'new_start'] = mult_df.loc[idx, 'new_start']
                            mult_df.loc[i, 'new_end'] = mult_df.loc[idx, 'new_end']
                            mult_df.loc[i, 'new_ID'] = mult_df.loc[idx, 'new_ID']
                            mult_df.loc[i, 'new_length'] = mult_df.loc[idx, 'new_length']
                            mult_df.loc[i, 'eval_gene'] = mult_df.loc[idx, 'eval_gene']
                            mult_df.loc[i, 'gene_chr'] = mult_df.loc[idx, 'ol_chr']
                            mult_df.loc[i, 'gene_start'] = mult_df.loc[idx, 'ol_start']
                            mult_df.loc[i, 'gene_end'] = mult_df.loc[idx, 'ol_end']
                            mult_df.loc[i, 'gene_ID'] = mult_df.loc[idx, 'ol_gene_ID'].split('/')[1]
                            mult_df.loc[i, 'gene_score'] = mult_df.loc[idx, 'ol_score']
                            mult_df.loc[i, 'gene_strand'] = mult_df.loc[idx, 'ol_strand']
                            mult_df.loc[i, 'gene_name'] = mult_df.loc[idx, 'ol_gene']
                            mult_df.loc[i, 'gene_feat_type'] = mult_df.loc[idx, 'ol_feat_type']
                            mult_df.loc[i, 'gene_length'] = mult_df.loc[idx, 'ol_gene_length']
                            if mult_df.loc[i, 'gene_strand'] == '+':
                                mult_df.loc[i, 'tes_start'] = mult_df.loc[i, 'gene_end'] - 1
                                mult_df.loc[i, 'tes_end'] = mult_df.loc[i, 'gene_end']
                            else:
                                mult_df.loc[i, 'tes_start'] = mult_df.loc[i, 'gene_start']
                                mult_df.loc[i, 'tes_end'] = mult_df.loc[i, 'gene_start'] + 1
                    continue
                # if no entries
                if pd.isna(mult_df.loc[idx, 'new_start']):
                    # change stranded column values
                    if mult_df.loc[idx, 'RT_strand'] == '+':
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[idx, 'RT_end']
                    else:
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[idx, 'RT_start']
                    # change non-stranded column values
                    mult_df.loc[idx, 'new_chr'] = mult_df.loc[idx, 'RT_chr']
                    mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                    mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    # make sure RT_length is a positive value (would not be if segment ends in middle of ds gene)
                    if mult_df.loc[idx, 'new_length'] < 0:
                        mult_df = mult_df.drop([idx])
                        continue
                    mult_df.loc[idx, 'num_geneOL'] = len(mult_df) - before_idx - 1
                    mult_df.loc[idx, 'RT_class_y'] = 'IId'
                    # overwrite original segment classes if small or not expressed
                    if (mult_df.loc[idx, 'ol_gene_length'] < 1000):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIa' + ',' + mult_df.loc[before_idx, 'RT_class_y']
                        continue
                    if (mult_df.loc[idx, 'expressed'] == 'no'):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIc' + ',' + mult_df.loc[before_idx, 'RT_class_y']
                        continue
                # could be IIf or IIIa; start of this could change, end shouldn't
                else:
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIf':
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                        else:
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])   
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIIa':
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                        else:
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start']) 
            # if bisect is in between, get index after and index before in list;
            # start is equal to start of nearest gene before that was evaled, end is equal to start of nearest ds
            # gene that was evaled unless it was IIf, in which case end does not change
            if (position != 0) & (position != len(index_list)):
                before_idx = index_list[position - 1]
                after_idx = index_list[position]
                # if idx is in the middle and equal to X, all subsequent rows should be changed to this unless 'IId'
                if mult_df.loc[idx, 'RT_class_y'] == 'X':
                    for i, r in mult_df.iloc[idx+1:len(mult_df)].iterrows():
                        if pd.isna(mult_df.loc[i, 'RT_class_y']):
                            mult_df.loc[i, 'new_chr'] = mult_df.loc[idx, 'new_chr']
                            mult_df.loc[i, 'new_start'] = mult_df.loc[idx, 'new_start']
                            mult_df.loc[i, 'new_end'] = mult_df.loc[idx, 'new_end']
                            mult_df.loc[i, 'new_ID'] = mult_df.loc[idx, 'new_ID']
                            mult_df.loc[i, 'new_length'] = mult_df.loc[idx, 'new_length']
                            mult_df.loc[i, 'eval_gene'] = mult_df.loc[idx, 'eval_gene']
                            mult_df.loc[i, 'gene_chr'] = mult_df.loc[idx, 'ol_chr']
                            mult_df.loc[i, 'gene_start'] = mult_df.loc[idx, 'ol_start']
                            mult_df.loc[i, 'gene_end'] = mult_df.loc[idx, 'ol_end']
                            mult_df.loc[i, 'gene_ID'] = mult_df.loc[idx, 'ol_gene_ID'].split('/')[1]
                            mult_df.loc[i, 'gene_score'] = mult_df.loc[idx, 'ol_score']
                            mult_df.loc[i, 'gene_strand'] = mult_df.loc[idx, 'ol_strand']
                            mult_df.loc[i, 'gene_name'] = mult_df.loc[idx, 'ol_gene']
                            mult_df.loc[i, 'gene_feat_type'] = mult_df.loc[idx, 'ol_feat_type']
                            mult_df.loc[i, 'gene_length'] = mult_df.loc[idx, 'ol_gene_length']
                            if mult_df.loc[i, 'gene_strand'] == '+':
                                mult_df.loc[i, 'tes_start'] = mult_df.loc[i, 'gene_end'] - 1
                                mult_df.loc[i, 'tes_end'] = mult_df.loc[i, 'gene_end']
                            else:
                                mult_df.loc[i, 'tes_start'] = mult_df.loc[i, 'gene_start']
                                mult_df.loc[i, 'tes_end'] = mult_df.loc[i, 'gene_start'] + 1
                    continue
                # if no entries
                if pd.isna(mult_df.loc[idx, 'new_start']):
                    if mult_df.loc[idx, 'RT_strand'] == '+':
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                    else:
                        mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                        mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                    # change non-stranded column values
                    mult_df.loc[idx, 'new_chr'] = mult_df.loc[idx, 'RT_chr']
                    mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                    mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    # make sure RT_length is a positive value (would not be if segment ends in middle of ds gene)
                    if mult_df.loc[idx, 'new_length'] < 0:
                        mult_df = mult_df.drop([idx])
                        continue
                    mult_df.loc[idx, 'num_geneOL'] = after_idx - before_idx -1
                    mult_df.loc[idx, 'RT_class_y'] = 'IId'
                    # overwrite original segment classes if small or not expressed
                    if (mult_df.loc[idx, 'ol_gene_length'] < 1000):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIa' + ',' + mult_df.loc[after_idx, 'RT_class_y']
                        continue
                    if (mult_df.loc[idx, 'expressed'] == 'no'):
                        mult_df.loc[idx, 'RT_class_x'] = 'IIc' + ',' + mult_df.loc[after_idx, 'RT_class_y']
                        continue
                # could be IIf or IIIa; either coordinate could change
                else: 
                    # could be IIf, this is incorrect, show warning and change to class IIIa
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIf':
                        warnings.warn("Exception: segment class IIf in middle position. Class changed to IIIa.")
                        mult_df.loc[idx, 'RT_class_y'] = 'IIIa'
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                        else:
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    if mult_df.loc[idx, 'RT_class_y'] == 'IIIa':
                        # change stranded column values
                        if mult_df.loc[idx, 'RT_strand'] == '+':
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[before_idx, 'ol_end']
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[after_idx, 'ol_start']
                        else:
                            mult_df.loc[idx, 'new_end'] = mult_df.loc[before_idx, 'ol_start']
                            mult_df.loc[idx, 'new_start'] = mult_df.loc[after_idx, 'ol_end']
                        # change non-stranded column values
                        mult_df.loc[idx, 'new_ID'] = mult_df.loc[idx, 'new_chr'] + ':' + str(int(mult_df.loc[idx, 'new_start'])) + '-' + str(int(mult_df.loc[idx, 'new_end']))
                        mult_df.loc[idx, 'new_length'] = int(mult_df.loc[idx, 'new_end'] - mult_df.loc[idx, 'new_start'])
                    
    return mult_df


def redefine_segment_coordinates(df, eval_df):
    '''
    Parameters
    ----------
    df: pandas dataframe containing reevaluated segment info (output from evaluate_segments())
    eval_df: pandas dataframe containing ambiguous segments from get_rt_eval_regions.py

    Description
    -------
    This function is the primary workflow function carrying out the reassignment of RT_segment coordinates after
    reevaluation (class IId). It first merges the reevaluated segment df (df) with the ambiguous segment df from
    get_rt_eval_regions.py (eval_df), and then reassigns RT_segments with single row entries in this df (no
    overlapping genes) and then segments with muli-row entries (1 or more overlapping genes). Output is a merged df
    with all redefined segment coordinates.
    
    Variables
    -------
    There are no variables that need to be changed in this function.
    '''
    
    # add ol_gene column
    eval_df['ol_gene'] = eval_df.ol_gene_ID.str.split('/', n=1, expand=True)[0]
    
    # change column name and merge df with eval df on original RT_ID, the ol_gene, and gene name
    df = df.rename(columns={'original_gene':'gene_name'})
    merged_df = pd.merge(eval_df, df, how='left', on=['RT_ID', 'ol_gene', 'gene_name', 'num_geneOL', 'expressed'])
    
    # fix the IIf class; should only apply for partial gene overlaps, if full overlap change to IIIa
    merged_df.loc[(merged_df['ol_gene_length'] == merged_df['ol_length']) & (merged_df['RT_class_y'] == 'IIf'), 'RT_class_y'] = 'IIIa'
    
    # rows in merged_df that could not be evaled are already in final_df, so don't need to worry about preserving
    # empty rows for single gene overlaps below
    
    # 1 gene overlaps can be set as RT_IDs directly, multi gene overlaps need to be parsed
    # split df to get segments that only overlapped 1 gene initially and those that overlapped multiple
    single = merged_df[merged_df['num_geneOL'] == 1]
    single = single.reset_index(drop=True)
    single['RT_class'] = single['RT_class_x'] + ',' + single['RT_class_y']
    
    fixed_single = single[single['new_start'].notnull()]
    fixed_single = fixed_single.reset_index(drop=True)
    fixed_single['original_RTseg'] = fixed_single['RT_ID']
    
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_chr'] = fixed_single.loc[fixed_single['RT_class'].str.contains('IId'),'new_chr']
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_start'] = fixed_single.loc[fixed_single['RT_class'].str.contains('IId'),'new_start']
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_end'] = fixed_single.loc[fixed_single['RT_class'].str.contains('IId'),'new_end']
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_ID'] = fixed_single.loc[fixed_single['RT_class'].str.contains('IId'),'new_ID']
    
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_score'] = '.'
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RT_length'] = fixed_single.loc[fixed_single['RT_class'].str.contains('IId'),'new_length']
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RTseg_counts'] = np.nan
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RTseg_density'] = np.nan
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'RTseg_rpkm'] = np.nan
    
    # will recompute on final df with bedtools, but doing this for now
    fixed_single.loc[fixed_single['RT_class'].str.contains('IId'), 'num_geneOL'] = 0
    
    # remove segment info (not accurate anymore) and new_ columns
    fixed_single_df = fixed_single[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class','RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm', 'original_RTseg',
                                    'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
                                    'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_coord', 'ol_feat_type', 'ol_gene_length', 'ol_length', 'num_geneOL', 'eval_gene', 'expressed']]
    
    # parse multi-gene overlaps
    multiple = merged_df[merged_df['num_geneOL'] > 1]
    multiple_set = set(multiple['RT_ID'])
    updated_df = pd.DataFrame()
    for seg in multiple_set:
        mult_df = multiple.loc[multiple['RT_ID'] == seg]
        # sort rows
        if all(mult_df['RT_strand'] == '+'):
            mult_df = mult_df.sort_values(by = ['ol_start'], ascending=True, ignore_index=True)
        else:
            mult_df = mult_df.sort_values(by = ['ol_end'], ascending=False, ignore_index=True)
        #check rows for class IId, and reassign RT_IDs where necessary
        if all(mult_df['RT_class_y'].isnull()):
            # add these to updated output, some not in final_df
            updated_df = pd.concat([updated_df, mult_df], ignore_index=True)
            continue
        if any(mult_df['RT_class_y'] == 'IIIa') and not any(mult_df['RT_class_y'] == 'IId'):
            # leave the segment as is, and add to updated_df
            updated_df = pd.concat([updated_df, mult_df], ignore_index=True)
            continue
        if any(mult_df['RT_class_y'] == 'IId'): # and not any(mult_df['RT_class_y'] == 'IIIa'):
            # change segment that were evaluated as unambiguous (class IId)
            # get indeces of any IId segments in mult_df
            index_list = mult_df.index[mult_df['RT_class_y'] == 'IId'].tolist()
            # add row below IId segment that represents the new segment created after the ol_gene
            for idx in index_list:
                # change stranded column values
                strand = mult_df.loc[idx, 'RT_strand']
                if strand == '+':
                    new_start = mult_df.loc[idx, 'ol_end']
                    new_end = mult_df.loc[idx, 'RT_end']
                else:
                    new_end = mult_df.loc[idx, 'ol_start']
                    new_start = mult_df.loc[idx, 'RT_start']
                # change non-stranded column values
                new_chr = mult_df.loc[idx, 'RT_chr']
                new_ID = str(new_chr) + ':' + str(int(new_start)) + '-' + str(int(new_end))
                new_length = int(new_end - new_start)
                # if new_length is negative, we don't want this new row, so continue
                if new_length < 0:
                    continue
                score = '.'
                # change RT_class_y to class I because this segment does not overlap anything
                RT_class_y = 'X'
                RT_feat_type = 'RT_segment'
                # add ol_gene cols because these will be gene cols
                ol_chr = mult_df.loc[idx, 'ol_chr']
                ol_start = mult_df.loc[idx, 'ol_start']
                ol_end = mult_df.loc[idx, 'ol_end']
                ol_gene_ID = mult_df.loc[idx, 'ol_gene_ID']
                ol_score = mult_df.loc[idx, 'ol_score']
                ol_strand = mult_df.loc[idx, 'ol_strand']
                ol_coord = mult_df.loc[idx, 'ol_coord']
                ol_feat_type = mult_df.loc[idx, 'ol_feat_type']
                ol_gene_length = mult_df.loc[idx, 'ol_gene_length']
                ol_gene = mult_df.loc[idx, 'ol_gene']
                # construct new row
                new_row = {'new_chr': new_chr,
                           'new_start': new_start,
                           'new_end': new_end,
                           'new_ID': new_ID,
                           'RT_score': score,
                           'RT_strand': strand,
                           'new_length': new_length,
                           'RT_class_y': RT_class_y,
                           'RT_feat_type': RT_feat_type,
                           'ol_chr': ol_chr,
                           'ol_start': ol_start,
                           'ol_end': ol_end,
                           'ol_gene_ID': ol_gene_ID,
                           'ol_score': ol_score,
                           'ol_strand': ol_strand,
                           'ol_coord': ol_coord,
                           'ol_feat_type': ol_feat_type,
                           'ol_gene_length': ol_gene_length,
                           'ol_gene': ol_gene,
                           'eval_gene': ol_gene}
                new_df = pd.DataFrame(new_row, index=[0])
                mult_df = pd.concat([mult_df.iloc[:idx+1], new_df, mult_df.iloc[idx+1:]], ignore_index=False, sort=False)
            # reset index after new row is added and get new index_list
            mult_df = mult_df.reset_index(drop=True)
            index_list = mult_df.index[mult_df['RT_class_y'] == 'IId'].tolist()
            # loop through rows of mult_df and if not in index_list, redefine segment based on other segments
            mult_df = update_segment_coordinates(mult_df, index_list)
            updated_df = pd.concat([updated_df, mult_df], ignore_index=True)
    
    # preserve information about original RT_segment for reference
    updated_df['original_RTseg'] = updated_df['RT_ID']
    updated_df['original_RTseg'].fillna(method='ffill', inplace=True)
    
    # fix RT_seg and gene info in rows that were changed (new_ cols become RT_ cols, ol_gene becomes gene)
    fixed_rows = []
    for idx, row in updated_df.iterrows():
        if pd.isna(row['new_start']):
            row['RT_class'] = row['RT_class_x']
            fixed_rows.append(row)
            continue
        else:
            if row['new_ID'] == row['RT_ID']:
                row['RT_class'] = row['RT_class_x']
                fixed_rows.append(row)
                continue
            else:
                # set RT_ columns to new columns
                row['RT_chr'] = row['new_chr']
                row['RT_start'] = row['new_start']
                row['RT_end'] = row['new_end']
                row['RT_ID'] = row['new_ID']
                if pd.isna(row['RT_class_x']):
                    row['RT_class'] = row['RT_class_y']
                if pd.isna(row['RT_class_y']):
                    row['RT_class'] = row['RT_class_x']
                if (pd.notna(row['RT_class_x'])) and (pd.notna(row['RT_class_y'])):
                    row['RT_class'] = str(row['RT_class_x']) + ',' + str(row['RT_class_y'])
                if 'IId' in row['RT_class']:
                    row['RT_score'] = '.'
                    row['RT_length'] = row['new_length']
                    row['RTseg_counts'] = np.nan
                    row['RTseg_density'] = np.nan
                    row['RTseg_rpkm'] = np.nan
                if (('IId' in row['RT_class']) and pd.notna(row['eval_gene'])) or pd.isna(row['gene_chr']):
                    if row['eval_gene'] != row['gene_name']:
                        # set ol_gene columns as gene columns
                        row['gene_chr'] = row['ol_chr']
                        row['gene_start'] = row['ol_start']
                        row['gene_end'] = row['ol_end']
                        row['gene_ID'] = row['ol_gene_ID'].split('/')[1]
                        row['gene_score'] = '.'
                        row['gene_strand'] = row['ol_strand']
                        row['gene_name'] = row['ol_gene']
                        row['gene_feat_type'] = row['ol_feat_type']
                        row['gene_length'] = row['new_length']
                        if row['gene_strand'] == '+':
                            row['tes_start'] = row['gene_end'] - 1
                            row['tes_end'] = row['gene_end']
                        else:
                            row['tes_start'] = row['gene_start']
                            row['tes_end'] = row['gene_start'] + 1
                fixed_rows.append(row)
    fixed_updated = pd.DataFrame(fixed_rows)
    
    # remove segment info (not accurate anymore) and new_ columns
    fixed_updated_df = fixed_updated[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class','RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm', 'original_RTseg',
                                      'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
                                      'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_coord', 'ol_feat_type', 'ol_gene_length', 'ol_length', 'num_geneOL', 'eval_gene', 'expressed']]
    
    # concat fixed_dfs and output this
    fixed_df = pd.concat([fixed_single_df, fixed_updated_df], ignore_index=True)
    fixed_df = fixed_updated_df.copy()
    # remove rows with genes who's tes does not match start of seg (mistakenly added somewhere)
    fixed_df = fixed_df[((fixed_df['RT_strand'] == '+') & (fixed_df['tes_end'] == fixed_df['RT_start'])) | ((fixed_df['RT_strand'] == '-') & (fixed_df['tes_start'] == fixed_df['RT_end']))]
    
    return fixed_df


def refine_classifications(final_df, fixed_df):
    
    '''
    Parameters
    ----------
    final_df: pandas dataframe containing all unambiguous segments from get_rt_eval_regions.py
    fixed_df: pandas dataframe containing fixed segment coordinates from reevaluated segments (output from
              redefine_segment_coordinates())
    
    Description
    -----------
    This function combines the fixed segments (fixed_df) with finalized segments from get_rt_eval_regions.py
    (final_df), removes redundancies within classifications, and assigns same classifications for any duplicate
    rows. The dataframe is then subsetted to contain only relevant columns and rows are deduplicated.
    
    Variables
    ---------
    
    '''
    
    # get list of original RT_IDs that were reevaluated in fixed_df to remove from final
    set_IDs = set(fixed_df['original_RTseg'])
    for seg in set_IDs:
        drop_list = final_df.index[final_df['RT_ID'] == seg].tolist()
        final_df = final_df.drop(labels=drop_list)
        
    final_df = final_df.reset_index(drop=True)
        
    # concatenate dfs, change column types, sort, and add original_RTseg column to empty rows (from final_df)
    df = pd.concat([final_df, fixed_df], ignore_index=True, sort=False)
    df[['RT_start', 'RT_end', 'tes_start', 'tes_end', 'gene_start', 'gene_end']] = df[['RT_start', 'RT_end', 'tes_start', 'tes_end', 'gene_start', 'gene_end']].astype('int64')
    df = df.sort_values(by = ['RT_chr', 'RT_start'], ignore_index=True)
    df.loc[pd.isna(df['original_RTseg']), 'original_RTseg'] = df.loc[pd.isna(df['original_RTseg']), 'RT_ID']
    
    
    class_row = []
    # remove redundant classes
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df)):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        # check for duplicates and remove them
        string = df.loc[idx, 'RT_class']
        new_classes = sorted(list(set(string.split(','))))
        # if contains IId and overlaps 0 genes, set RT_class = I,IId
        # if contains IId and IIb, remove IIb (cannot partially overlap anything)
        if 'IId' in new_classes:
            if row['num_geneOL'] == 0:
                row['RT_class'] = 'I,IId'
                class_row.append(row)
                continue
            if 'IIb' in new_classes:
                new_classes.remove('IIb')
        # if contains IIb and IIe/IIf, only needs IIe/IIf
        if any(x == 'IIb' for x in new_classes) and any((x == 'IIe' or x == 'IIf') for x in new_classes):
            new_classes.remove('IIb')
        # join back together and set as RT_class
        new_string = ','.join(sorted(new_classes, key=len))
        row['RT_class'] = new_string
        class_row.append(row)
    class_df = pd.DataFrame(class_row)
    
    # make sure that all new RT_IDs have the same classes in all rows
    new_ids = set(class_df['RT_ID'])
    for seg in new_ids:
        df_sub = class_df[class_df['RT_ID'] == seg]
        df_sub = df_sub.reset_index(drop=True)
        if len(df_sub) == 1:
            continue
        if all(df_sub['RT_class'] == df_sub['RT_class'][0]):
            continue
        else:
            # get list of all RT_class values, split str, and flatten list
            string_list = list(df_sub['RT_class'])
            class_list = [x.split(',') for x in string_list]
            class_list = [i for j in class_list for i in j]
            # get set of uniqe values as list and join (sort alaphabetically then by length to get correct order)
            class_set = sorted(list(set(class_list)))
            RT_class = ','.join(sorted(class_set, key=len))
            # set new class string as RT_class for all rows
            class_df.loc[class_df['RT_ID'] == seg, 'RT_class'] = RT_class
            
    # subset df to just contain RT_seg and gene columns and remove duplicate rows
    output_df = class_df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class','RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm', 'original_RTseg',
                          'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm'
                          ]].copy()
    
    output_df['RT_length'] = output_df['RT_end'] - output_df['RT_start']
    output_df = output_df.round(decimals={'RTseg_counts': 5, 'RTseg_density': 5, 'RTseg_rpkm': 5, 'gene_counts':5, 'gene_density':5, 'gene_rpkm':5})
    
    output_df_dedup = output_df.drop_duplicates(ignore_index=True)
    
    return output_df_dedup


def main(start_file, end_file, gene_file, final_file, eval_file):
    
    '''
    Parameters
    ----------
    start_file: [arg1] 
    end_file: [arg2]
    gene_file: [arg3]
    final_file: [arg4]
    
    Description
    -----------
    This function parses ambiguous segments to determine which can be evaluated as unambiguous and then refines 
    classifications of all segments. Classifications are required for evaluation, but will be reassigned in the 
    final step.

    
    Variables
    ---------
    No variables to change in this main function.
    '''
    
    start_df = load_files(start_file)
    end_df = load_files(end_file)
    gene_df = load_files(gene_file)
    final_df = load_files(final_file)
    eval_df = load_files(eval_file)
    
    # parse eval segments to determine which can be evaluated as unambiguous and classify all segments
    evaluated_seg = evaluate_segments(start_df, end_df, gene_df)
    
    # merge evaluated segments back with the rest of their info from the original df and modify segment coordinates
    fixed_df = redefine_segment_coordinates(evaluated_seg, eval_df)
    
    # refine previous classifications
    output_df = refine_classifications(final_df, fixed_df)
    
    # output final df
    output_df = output_df.sort_values(by = ['RT_chr', 'RT_start'])
    output_df = output_df.rename(columns={"RT_chr":"#RT_chr"})
    
    # output file (output should be: $SAMPLE.RTsegments.TESintersect.segmentRPKM.fixed.bed)
    output_string = os.path.splitext(final_file)[0]
    fixed_output = output_string.rsplit(sep='.', maxsplit=1)[0]+".fixed.bed"
    delete_if_exists(fixed_output)
    output_df.to_csv(fixed_output, sep='\t', header=True, index=False)
    


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('eval_start_file',
                        help='Path to eval file containing RT start counts (expects .segmentRPKM.startsRPKM.bed file).')
    parser.add_argument('eval_end_file',
                        help='Path to eval file containing RT end counts (expects .segmentRPKM.endsRPKM.bed file).')
    parser.add_argument('eval_gene_file',
                        help='Path to eval file containing downstream gene counts (expects .segmentRPKM.genesRPKM.bed file).')
    parser.add_argument('noeval_file',
                        help='Path to file containing unambiguous, previously classified segments (expects .segmentRPKM.noeval.bed).')
    parser.add_argument('eval_file',
                        help='Path to file containing ambiguous, previously classified segments (expects .segmentRPKM.eval.bed).')

    
    args = parser.parse_args()

    main(args.eval_start_file, args.eval_end_file, args.eval_gene_file, args.noeval_file, args.eval_file)
    
    # start_file = None
    # end_file = None
    # gene_file = None
    # final_file = None
    # eval_file = None
    
    # main(start_file, end_file, gene_file, final_file, eval_file)