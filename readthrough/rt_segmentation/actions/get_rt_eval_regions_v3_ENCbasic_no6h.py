# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 17:11:28 2022

Description
-----------
This script defines the following functions which are used to... :
    
    load_seg_files: loads in segment files needed for this script
    load_exp_files: loads in expression files needed for this script
    delete_if_exists: deletes file if already exists in working directory
    format_for_final_df: reformats dfs for the final_df output
    is_overlap: checks if there is any genes overlapped by an RT_segment
    classify_gene_overlap: classifies the basic gene overlaps (types IIa, IIb, or IIc)
    sort_partial_genes: identifies ambiguous segments where downstream gene is overlapping with another gene
    check_6h_expression: determines if downstream genes are expressed in 6h data
    find_short_overlaps: classifies gene overlaps that are less than 1kb
    classify_segments: master script that carries out segment classifications
    redefine_upstream_gene: redefines upstream gene for eval (relevant for segments overlapping more than 1 gene)
    isNaN: assesses nan values when they are not numerical (np.nan will not work)
    get_eval_regions: master script that defines upstream genes and obtains regions for eval

Version Notes
-------------
v2: modified methods and removed excess code
v3: removed excess code and changed short gene overlap threshold to 2kb

@author: abmcs
"""


import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import warnings


def load_seg_files(filepath):
    
    '''
    loads in pandas dataframes of each of the segment files required for this step, adds a header, adds 3 new
    columns: 1.) containing the coordinate ID for the overlapping genes found ('ol_coord' --> chr:start-end), 
    2.) containing the gene_length (bp) of the overlapping genes found, 3.) containing the number of gene overlaps
    for each RT_segment, and restructures the column orgnaization of the df for use in downstream steps (columns
    removed and reorganized)
    '''
    
    df = pd.read_csv(filepath, sep='\t', names=None, engine='python')
    
    # adding this if statement because of count_bed weirdness --> sometimes it erases end col of some files when
    # ouputting, this helps avoid errors from extra col added to circumvent that
    if df.iloc[0,42] == 'X':
        header = ['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_length',
              'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
              'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'segment_end_score', 'num_fused_seg', 'intersect_bp',
              'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_length', 'X',
              'seglen', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm']
    
        df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    else:
        header = ['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_length',
              'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
              'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'segment_end_score', 'num_fused_seg', 'intersect_bp',
              'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_length',
              'seglen', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm']
    
        df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    
    # create coordinate ID column for overlapping gene
    df['ol_coord'] = df['ol_chr'].astype(str) + ":" + df['ol_start'].astype(str) + "-" + df['ol_end'].astype(str)
    
    # create column with overlapping gene length
    df['ol_gene_length'] = df['ol_end'].astype(int) - df['ol_start'].astype(int)
    
    # create column with counts of gene overlaps per segment (segs with no overlap will say 1, fixed later)
    counts = df['RT_ID'].value_counts().to_frame()
    counts.reset_index(inplace=True)
    counts.rename(columns={'RT_ID':'num_geneOL', 'index':'RT_ID'}, inplace=True)
    df = pd.merge(df, counts, on=['RT_ID'], how='outer')
    
    # reformat df
    formatted_df = df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm',
                       'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
                       'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'segment_end_score', 'num_fused_seg',
                       'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_coord', 'ol_feat_type', 'ol_gene_length', 'ol_length', 'num_geneOL']]

    return formatted_df


def load_exp_files(filepath, sample, assay=None):
    
    '''
    loads in pandas dataframes of the file containing gene expression values required for this step,
    creates two new columns with gene name (gene_name) and ensembl id (gene_ID) separated, and subsets df to
    contain only rpkm values for sample
    '''
    
    df = pd.read_csv(filepath, sep='\t', engine='python')
    df = df.rename(columns={'geneName/geneID': 'gene'})
    
    if assay == '6h':
        sample_6h = sample.replace('0h', '6h')
        rpkm_col = 'exon_RPKM_' + sample_6h
    else:
        rpkm_col = 'gene_RPKM_' + sample
    
    df[['gene_name', 'gene_ID']] = df['gene'].str.rsplit('/', n=1, expand=True)
    
    rpkm_sub = df[['gene_ID', 'gene_name', rpkm_col]]
    
    out_df = rpkm_sub.rename(columns={rpkm_col:'rpkm'})
        
    return out_df


def delete_if_exists(filepath):
    
    '''
    removes existing output file if filename matches one in directory
    '''
    
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    
    return False


def format_for_final_df(df):
    
    '''
    removes overlapping gene columns ('ol_') and rearranges column order for final df output
    '''
    
    dropped_df = df.drop(list(df.filter(regex = 'ol_')), axis = 1)
    
    formatted_df = dropped_df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm',
                               'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
                               'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'segment_end_score', 'num_fused_seg', 'num_geneOL']]
    
    return formatted_df


def is_overlap(df):
    
    '''
    checks each row of input/working df to check if there is an overlapping gene recorded for each RT_segment and
    adds these to two new lists to output as dfs with either segments that overlap at least 1 gene (yes_df) or
    segments that overlap no genes (no_df)
    '''
    
    yes_overlap = []
    no_overlap = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Check for gene overlaps"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        if row['ol_chr'] != '.':
            yes_overlap.append(row)
        else:
            row['num_geneOL'] = 0
            no_overlap.append(row)
    
    yes_df = pd.DataFrame(yes_overlap)
    
    no_df = pd.DataFrame(no_overlap)
    no_df['RT_class'] = 'I'
    
    return yes_df, no_df


def classify_gene_overlap(df):
    
    '''
    Parameters
    ----------
    df: pandas df containing RT_segments that overlap at least one gene
    
    Description
    -----------
    This function evaluates each RT_segments gene overlap status to determine the types of genes overlapped by the
    segment, and assigns them the cooresponding classification.
    
    Variables
    ---------
    No variables require changing, unless classes are to be renamed
    '''

    df = df.reset_index(drop=True)
    
    # add RT_class column if needed
    if 'RT_class' not in df:
        df['RT_class'] = np.nan
    
    # short
    short_ol_df = df[df['ol_gene_length'] < 1000]
    short_ol = list(short_ol_df.RT_ID.unique())
    for seg in short_ol:
        df.loc[df['RT_ID'] == seg, 'RT_class'] = 'IIa'
    
    # long
    full_ol_df = df[df['ol_gene_length'] >= 1000]
    full_ol_df = full_ol_df[full_ol_df['ol_gene_length'] == full_ol_df['ol_length']]
    full_ol = list(full_ol_df.RT_ID.unique())
    for seg in full_ol:
        if all(df.loc[df['RT_ID'] == seg, 'RT_class'].isnull()):
            df.loc[df['RT_ID'] == seg, 'RT_class'] = 'IIc' # could also be 'IId,III' technically, but will add these later
        else:
            df.loc[df['RT_ID'] == seg, 'RT_class'] = df.loc[df['RT_ID'] == seg, 'RT_class']+',IIc' # could be 'IId,III'
    
    # partial
    part_df = df.loc[df['ol_gene_length'] >= 1000]
    part_df = part_df[part_df['ol_gene_length'] != part_df['ol_length']]
    part_ol = list(part_df.RT_ID.unique())
    for seg in part_ol:
        if all(df.loc[df['RT_ID'] == seg, 'RT_class'].isnull()):
            df.loc[df['RT_ID'] == seg, 'RT_class'] = 'IIb' # could also be 'IIe,III' technically, but will add these later
        else:
            df.loc[df['RT_ID'] == seg, 'RT_class'] = df.loc[df['RT_ID'] == seg, 'RT_class']+',IIb' # could be 'IIe,III'
    
    return df


def sort_partial_genes(df):
    
    '''
    Parameters
    ----------
    df: pandas df containing RT_segments that overlap at least one gene that is greater than 1kb in length
    
    Description
    -----------
    This function classifies RT_segments that have a partial overlap with a gene due to the gene of origin and the
    "downstream" gene overlapping to some degree. These RT_segments are not able to be evaluated and will remain
    ambiguous.
    
    Variables
    ---------
    No variables require changing, unless classes are to be renamed
    '''
    
    df = df.reset_index(drop=True)
    
    ambiguous_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Find ambiguous segments"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        # verify that start of ol_gene is downstream of end of RT gene, if not add to ambiguous_rows and remove from df
        if row['RT_strand'] == '+':
            if  row['gene_end'] > row['ol_start']:
                ambiguous_rows.append(row)
                #df.drop([idx], inplace = True)
        else:
            if row['gene_start'] < row['ol_end']:
                ambiguous_rows.append(row)
                #df.drop([idx], inplace = True)
                
    ambiguous_df = pd.DataFrame(ambiguous_rows)
    ambiguous_df = ambiguous_df.reset_index(drop=True)
    
    # classify ambiguous segments
    # for row_ind, row_vals in tqdm(enumerate(ambiguous_df.iterrows()), total=len(ambiguous_df), desc="Classify ambiguous segments"):
    for row_ind, row in ambiguous_df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        if row['RT_strand'] == '+':
            # check if gene is completely overlapped by ol_gene
            if (row['gene_start'] >= row['ol_start']) and (row['gene_end'] <= row['ol_end']):
                ambiguous_df.loc[idx, 'RT_class'] = 'IIIb'
            # check if gene is partially overlapped by ol_gene
            if (row['gene_start'] < row['ol_start']) and (row['gene_end'] > row['ol_start']):
                ambiguous_df.loc[idx, 'RT_class'] = 'IIIc'
        else:
            # check if gene is completely overlapped by ol_gene
            if (row['gene_start'] >= row['ol_start']) and (row['gene_end'] <= row['ol_end']):
                ambiguous_df.loc[idx, 'RT_class'] = 'IIIb'
            # check if gene is partially overlapped by ol_gene
            if (row['gene_end'] > row['ol_end']) and (row['gene_start'] < row['ol_end']):
                ambiguous_df.loc[idx, 'RT_class'] = 'IIIc'             
   
    # populate original df with these classes
    ambig_ol = list(ambiguous_df.RT_ID.unique())
    for seg in ambig_ol:
        new_class = ambiguous_df.loc[ambiguous_df['RT_ID'] == seg, 'RT_class']
        if all(df.loc[df['RT_ID'] == seg, 'RT_class'].isnull()):
            if len(new_class) == 1:
                (df.loc[df['RT_ID'] == seg, 'RT_class']) = new_class.item()
            else:
                new_classes = set(new_class)
                new_classes = ','.join(new_classes)
                (df.loc[df['RT_ID'] == seg, 'RT_class']) = new_classes
        else:
            if len(new_class) == 1:
                (df.loc[df['RT_ID'] == seg, 'RT_class']) = (df.loc[df['RT_ID'] == seg, 'RT_class'])+','+new_class.item()
            else:
                new_classes = set(new_class)
                new_classes = ','.join(new_classes)
                (df.loc[df['RT_ID'] == seg, 'RT_class']) = (df.loc[df['RT_ID'] == seg, 'RT_class'])+','+new_classes
    
    # subset final df into ambiguous or not
    ambiguous_final = df[df['RT_class'].str.contains('IIIb|IIIc')]
    df_final = df[~df['RT_class'].str.contains('IIIb|IIIc')]
    
    return df_final, ambiguous_final


def check_6h_expression(df, df_exon_exp):
    
    '''
    Parameters
    ----------
    df: pandas df containing RT_segments that overlap at least one gene that is greater than 1kb in length and
        where gene overlaps are not ambiguous (output from sort_partial_genes)
    df_exon_exp: pandas df containing BruChase exon RPKM values for evaluating gene expression at 6h timepoint
    
    Description
    -----------
    This function checks the expression of each overlapping gene at 6h and adds an "expressed" column to df. If 
    the gene is not expressed, the RT_segment is designated as such. If the gene is expressed (RPKM > 0.5), the 
    segment will need to be evaluated further.
    
    Variables
    ---------
    6h gene expression cutoff can be modified (currently >0.5) --> modified to 0.25 becuase this script uses 0h gene rpkm
    '''
    
    df = df.reset_index(drop=True)
    
    # add expressed column
    df['expressed'] = np.nan
    
    expressed_df_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Check for expression in 6h data"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        gene_name = df.loc[idx, 'ol_gene_ID'].split('/')[0]
        exp = df_exon_exp.loc[df_exon_exp['gene_name'] == gene_name, 'rpkm']
        # make sure there is a gene exp value, if not remove row
        if exp.empty:
            row['expressed'] = 'missing'
            continue
        else:
            exp = exp.values[0]
        
        if exp > 0.25:
            row['expressed'] = 'yes'
            expressed_df_rows.append(row)
        else:
            row['expressed'] = 'no'
            expressed_df_rows.append(row)
    
    expression_df = pd.DataFrame(expressed_df_rows)
    
    # create missing_df for rows with no expression value
    missing_df = expression_df[expression_df['expressed'] == 'missing']
    expression_df = expression_df[expression_df['expressed'] != 'missing']
    
    # create not_expressed_df (segment where all genes are not expressed)
    not_expressed_df = pd.DataFrame()
    expressed_df = pd.DataFrame()
    id_list = list(expression_df.RT_ID.unique())
    for seg in id_list:
        new_df = expression_df[expression_df['RT_ID'] == seg]
        if all(new_df['expressed'] == 'no'):
            not_expressed_df = pd.concat([not_expressed_df, new_df], ignore_index=True)
        else:
            expressed_df = pd.concat([expressed_df, new_df], ignore_index=True)
    
    return expressed_df, not_expressed_df, missing_df


def find_short_overlaps(df):
    
    '''
    Parameters
    ----------
    df: pandas df containing RT_segments that overlap at least one gene that is greater than 1kb in length, where
        gene overlaps are not ambiguous, and where gene is expressed at 6h (output from check_6h_expression)
    
    Description
    -----------
    This function checks the degree of overlap of each RT_segment with the expressed downstream gene. If the
    overlap with the downstream gene is less than 1kb, this segment is classified as such. If the
    overlap is greater than 1kb, then the segment will require further evaluation.
    
    Variables
    ---------
    Overlap distance with downstream gene can be modified (currently >2000)
    '''
    
    df = df.reset_index(drop=True)
    
    #short_ol_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Find short overlaps"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        #idx = row_ind
        # do not evaluate a gene that is less than 1kb or that is overlapped entirely by RT_seg
        if row['ol_gene_length'] < 2000:
            continue
        if row['ol_gene_length'] == row['ol_length']:
            continue
        # check if end of segment is within 2kb of TSS, if yes append to short_ol_rows
        if row['ol_length'] < 2000:
            seg = row['RT_ID']
            seg_rows = df.loc[df['RT_ID'] == seg]
            if len(seg_rows) == 1:
                df.loc[df['RT_ID'] == seg, 'RT_class'] = 'IIe'
            else:
                df.loc[df['RT_ID'] == seg, 'RT_class'] = df.loc[df['RT_ID'] == seg, 'RT_class'] + ',IIe'
    
    # get segments where this is the only overlap and remove from eval_df output
    short_ol_df = df[df['RT_class'] == 'IIe']
    eval_df = df[df['RT_class'] != 'IIe']
    
    return eval_df, short_ol_df


def classify_segments(df, df_exon_exp):
    
    '''
    Parameters
    ----------
    df: pandas df of primary working file (.segmentsRPKM.bed) that contains all RT_segments, all segment gene
        overlaps, and all segment RPKM values
    df_exon_exp: pandas df containing BruChase exon RPKM values for evaluating gene expression at 6h timepoint
    
    Description
    -----------
    This function runs the necessary functions to classify RT_segments based on gene overlaps and determine if
    further evaluation is needed. The outputs include a pandas dataframe with segments that are still ambiguous
    and require downstream evaluation (eval_segments) and a pandas dataframe containing all successfully classified
    segments (final_df). Also outputted is a dataframe containing segments with overlapping genes that were not
    present in the 6h gene expression df (missing_df).
    
    Variables
    ---------
    6h expression cutoff can be modified in check_6h_expression(), and short overlap cutoff can be modified in
    find_short_overlaps() if desired
    '''
    
    df = df.reset_index(drop=True)
    
    # First, find and designate all classes in df
    # get rows with no gene overlaps, and add no_overlap to final_df
    continue_eval, no_overlap = is_overlap(df)
    final_df = format_for_final_df(no_overlap) # no need to dedup no_overlap
    
    eval_df = classify_gene_overlap(continue_eval)
    
    # Second, find segments that only overlap short genes, add to final_list
    only_short = eval_df[eval_df['RT_class'] == 'IIa']
    only_short = format_for_final_df(only_short)
    only_short_dedup = only_short.drop_duplicates(subset=['RT_ID'], ignore_index=True)
    final_df = pd.concat([final_df, only_short_dedup], ignore_index=True)
    
    # Third, focus on long genes and determine which should be evaluated further
    # remove any overlaps of partial genes where end of gene does not fall in segment
    eval_df2 = eval_df[eval_df['RT_class'] != 'IIa']
    
    eval_df3, ambiguous_df = sort_partial_genes(eval_df2)
    
    # add ambiguous_df to final df, because these probably cannot be parsed further
    ambiguous_df = format_for_final_df(ambiguous_df)
    ambiguous_df_dedup = ambiguous_df.drop_duplicates(subset=['RT_ID'], ignore_index=True)
    final_df = pd.concat([final_df, ambiguous_df_dedup], ignore_index=True)
    
    # Check ol genes for expression at 6h; if rpkm < 0.5, add to final_list
    expressed, not_expressed, missing_exp = check_6h_expression(eval_df3, df_exon_exp)
    not_expressed = format_for_final_df(not_expressed)
    not_expressed_dedup = not_expressed.drop_duplicates(subset=['RT_ID'], ignore_index=True)
    final_df = pd.concat([final_df, not_expressed_dedup], ignore_index=True)
    
    # For expressed genes, check if overlap between segment and gene is less than 1kb, if yes classify as such
    eval_segments, small_overlaps = find_short_overlaps(expressed)
    small_overlaps = format_for_final_df(small_overlaps)
    small_overlaps_dedup = small_overlaps.drop_duplicates(subset=['RT_ID'], ignore_index=True) # should not be duplicates because these do not overlap more than one gene
    final_df = pd.concat([final_df, small_overlaps_dedup], ignore_index=True)

    return eval_segments, final_df, missing_exp


def redefine_upstream_gene(dup_df, fixed_rows):
    
    if all(dup_df['RT_strand'] == '+'):
        dup_df = dup_df.sort_values(by = ['ol_start'], ascending=True, ignore_index=True)
        for idx, row in dup_df.iterrows():
            if idx != 0:
                # set previous row based on index
                prev_row = dup_df.iloc[idx-1]
                if prev_row['ol_end'] < row['ol_start']:
                    gene_name = prev_row['ol_gene_ID'].split('/')[0]
                    # replace eval gene info with that of the previous row
                    row['eval_chr'] = prev_row['ol_chr']
                    row['eval_start'] = int(prev_row['ol_start'])
                    row['eval_end'] = int(prev_row['ol_end'])
                    row['eval_gene'] = gene_name
                    row['eval_score'] = prev_row['ol_score']
                    row['eval_strand'] = prev_row['ol_strand']
                    row['eval_feat_type'] = prev_row['ol_feat_type']
                else:
                    row['eval_chr'] = np.nan
                    row['eval_start'] = np.nan
                    row['eval_end'] = np.nan
                    row['eval_gene'] = np.nan
                    row['eval_score'] = np.nan
                    row['eval_strand'] = np.nan
                    row['eval_feat_type'] = np.nan
                # add to fixed_rows
                fixed_rows.append(row)
            else:
                # if index is zero, eval gene does not change
                fixed_rows.append(row)
    else:
        dup_df = dup_df.sort_values(by = ['ol_end'], ascending=False, ignore_index=True)
        for idx, row in dup_df.iterrows():
            if idx != 0:
                # set previous row based on index
                prev_row = dup_df.iloc[idx-1]
                if prev_row['ol_start'] > row['ol_end']:
                    gene_name = prev_row['ol_gene_ID'].split('/')[0]
                    # replace eval gene info with that of the previous row
                    row['eval_chr'] = prev_row['ol_chr']
                    row['eval_start'] = int(prev_row['ol_start'])
                    row['eval_end'] = int(prev_row['ol_end'])
                    row['eval_gene'] = gene_name
                    row['eval_score'] = prev_row['ol_score']
                    row['eval_strand'] = prev_row['ol_strand']
                    row['eval_feat_type'] = prev_row['ol_feat_type']
                else:
                    row['eval_chr'] = np.nan
                    row['eval_start'] = np.nan
                    row['eval_end'] = np.nan
                    row['eval_gene'] = np.nan
                    row['eval_score'] = np.nan
                    row['eval_strand'] = np.nan
                    row['eval_feat_type'] = np.nan
                # add to fixed_rows
                fixed_rows.append(row)
            else:
                # if index is zero, eval gene does not change
                fixed_rows.append(row)
    
    return fixed_rows


def isNaN(val):
    return val != val


def get_eval_regions(df):
    
    '''
    Parameters
    ----------
    df: pandas df containing RT_segments that are not able to be classified and need to be more rigorously
    evaluated downstream (output from classify_segments)
    
    Description
    -----------
    This function creates three new dfs (bed format) containing the coordinates of regions (input for count_bed)
    for downstream evaluation:
        1. start_df: beginning of RT_segment (RT_seg start plus 1kb)
        2. end_df: end of RT_segment (gene start minus 1kb)
        3. gene_df: beginngin of downstream gene (gene start plus 1kb)
    
    Variables
    ---------
    Evaluation region size can be modified (currently 1000, thus distance between genes is >2000)
    '''
    
    df = df.reset_index(drop=True)
    
    # remove rows with small genes (ignored in downstream eval) and with small overlaps and no expression
    df_subset = df[df['ol_gene_length'] > 1000]
    df_subset = df_subset[df_subset['ol_length'] > 1000]
    df_subset = df_subset[df_subset['expressed'] == 'yes']
    
    # add eval columns to gather gene information for the upstream gene
    df_subset['eval_chr'] = df_subset['gene_chr']
    df_subset['eval_start'] = df_subset['gene_start']
    df_subset['eval_end'] = df_subset['gene_end']
    df_subset['eval_gene'] = df_subset['gene_name']
    df_subset['eval_score'] = df_subset['gene_score']
    df_subset['eval_strand'] = df_subset['gene_strand']
    df_subset['eval_feat_type'] = df_subset['gene_feat_type']
    
    # separate out duplicates and reassign the upstream gene to the nearest expressed gene
    duplicate_df = df_subset[df_subset['RT_ID'].duplicated(keep=False)]
    duplicates = duplicate_df.RT_ID.unique()
    
    fixed_rows = []
    for seg in duplicates:
        dup_df = duplicate_df[duplicate_df['RT_ID'] == seg]
        fixed_rows = redefine_upstream_gene(dup_df, fixed_rows)
    
    fixed_df = pd.DataFrame(fixed_rows)
    fixed_df = fixed_df.reset_index(drop=True)
    
    # fix the stupid scientific notation
    fdf_start = fixed_df['eval_start'].apply(lambda x: '%.0f' % x)
    fixed_df['eval_start'] = fdf_start
    fdf_end = fixed_df['eval_end'].apply(lambda x: '%.0f' % x)
    fixed_df['eval_end'] = fdf_end
    
    unique_df = df_subset.drop_duplicates(subset=['RT_ID'], keep=False)
    modified_df = pd.concat([unique_df, fixed_df], ignore_index=True)
    
    # separate segments based on distance between genes and create eval_df
    start_rows = [] # RT_start
    end_rows = [] # RT_end
    gene_rows = [] # downstream_gene_start
    ambiguous_rows = []
    # for row_ind, row_vals in tqdm(enumerate(modified_df.iterrows()), total=len(modified_df)):
    for row_ind, row in modified_df.iterrows():
        # row = row_vals[1]
        # verify that strands of segment and ol_genes match
        if row['RT_strand'] != row['ol_strand']:
            error_message = "ERROR [RT_ID: "+str(row["RT_ID"])+"]: Segment and Gene are not on the same strand"
            sys.exit(error_message)
        
        if isNaN(row['eval_chr']):
            ambiguous_rows.append(row)
            continue
    
        if row['eval_strand'] == '+':
            eval_end = int(row['eval_end'])
            ol_start = row['ol_start']
            gene_dist = ol_start - eval_end
            if gene_dist > 2000:
                RTst_start = eval_end
                RTst_end = eval_end + 1000
                RTen_start = ol_start - 1000
                RTen_end = ol_start
                DSG_start = ol_start
                DSG_end = ol_start + 1000
                ol_gene = row['ol_gene_ID'].split('/')[0]
                
                start_rows.append([row['eval_chr'], RTst_start, RTst_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
                end_rows.append([row['eval_chr'], RTen_start, RTen_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
                gene_rows.append([row['eval_chr'], DSG_start, DSG_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
            else:
                ambiguous_rows.append(row)
        else:
            eval_start = int(row['eval_start'])
            ol_end = row['ol_end']
            gene_dist = eval_start - ol_end
            if gene_dist > 2000:
                RTst_start = eval_start - 1000
                RTst_end = eval_start
                RTen_start = ol_end
                RTen_end = ol_end + 1000
                DSG_start = ol_end - 1000
                DSG_end = ol_end
                ol_gene = row['ol_gene_ID'].split('/')[0]
                
                start_rows.append([row['eval_chr'], RTst_start, RTst_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
                end_rows.append([row['eval_chr'], RTen_start, RTen_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
                gene_rows.append([row['eval_chr'], DSG_start, DSG_end, row['RT_ID'], row['eval_score'], row['eval_strand'], row['gene_name'], row['eval_gene'], ol_gene, gene_dist, row['expressed'], row['num_geneOL']])
            else:
                ambiguous_rows.append(row)
            
    start_df = pd.DataFrame(start_rows, columns = ['chr', 'start', 'end', 'RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL'])
    end_df = pd.DataFrame(end_rows, columns = ['chr', 'start', 'end', 'RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL'])
    gene_df = pd.DataFrame(gene_rows, columns = ['chr', 'start', 'end', 'RT_ID', 'score', 'strand', 'original_gene', 'eval_gene', 'ol_gene', 'dist', 'expressed', 'num_geneOL'])
    ambiguous_df = pd.DataFrame(ambiguous_rows)
    
    # add ambiguous RT_class to genes in ambiguous df (too close together to evaluate)
    ambiguous_df['RT_class'] = ambiguous_df['RT_class'].astype(str) + ',IIIa'
    
    # set ambiguous_df as final_df
    final_df = ambiguous_df
    
    # get small genes, small_dist, and not expressed segs, and add to a new df
    small_genes = df[df['ol_gene_length'] < 1000]
    small_dist = df[df['ol_length'] < 1000]
    not_expressed = df[df['expressed'] == 'no']
    no_eval_df = pd.concat([small_genes, small_dist], ignore_index=True)
    no_eval_df = pd.concat([no_eval_df, not_expressed], ignore_index=True)
    
    # TEST: check overlaps of no_eval_df and subsetted df
    # (there should be no RT_seg in no_eval that are not found in the subsetted eval df)
    uniq_ne = set(no_eval_df['RT_ID'])
    uniq_sub = set(df_subset['RT_ID'])
    diff =  uniq_ne - uniq_sub
    if len(diff) != 0:
        # add no_eval segments that are not in the eval subset to final_df along with ambiguous_df
        # add rest to other and confirm that they are all present in the subset, else warning
        other_df = pd.DataFrame()
        for seg in uniq_ne:
            #seg = "chr1:173858808-173868500"
            rows = df.loc[df['RT_ID'] == seg]
            if all(map(lambda x, y: x or y, (rows['ol_length'].values < 1000), (rows['expressed'].values == 'no'))):
                final_df = pd.concat([final_df, rows], ignore_index=True)
            else:
                other_df = pd.concat([other_df, rows], ignore_index=True)
    
    # check again and confirm these other segments are represented in df_sub
    uniq_oth = set(other_df['RT_ID'])
    diff2 = uniq_oth - uniq_sub
    if len(diff2) != 0:
        warnings.warn("Segments are missing in eval. Check output of get_eval_regions().")
        
    return start_df, end_df, gene_df, final_df


def main(working_file, sample, exon_exp_file):
    
    '''
    Parameters
    ----------
    working_file: [arg1] path to RT_seg .bed file with gene overlap and count_bed info (expects .segmentRPKM.bed file)
    exon_exp_file: [arg2] path to matrix containing 6h exon RPKM values
    
    Description
    -----------
    Loads in RTsegment file (must contain TES info, proper segment gene overlaps, and segment count/rpkm values),
    classifies RTsegments, filters the segments to get those that need to be evaluated further, and produces new
    files containing the coordinates of the eval regions that will be inputs for count_bed.
    
    Variables
    ---------
    No variables to change in this main function.
    '''
    
    if working_file is None:
        working_file = "C:/Users/abmcs/Desktop/shared_biosamples/read-through/fused_dec/bin_250_old/HCT1160h4001a.RTsegments.TESintersect.segmentRPKM.bed"
    # # build in workaround for no 6h data
    # if exon_exp_file is None:
        # exon_exp_file = gene_exp_file
        # assay = None
    # else:
        # assay = '6h'
    
    
    df_seg = load_seg_files(working_file)
    df_exon_exp = load_exp_files(exon_exp_file, sample, assay=None)
    
    # take RTseg count cutoff
    df_cutoff = df_seg[df_seg['RTseg_counts'] > 10]
    
    # classify segments and get subset for eval
    # get segments ol 0 genes, fix formatting, and output as final_df
    eval_segments, final_df, missing_genes = classify_segments(df_cutoff, df_exon_exp)
    
    # get eval regions
    eval_start, eval_end, eval_gene, no_eval = get_eval_regions(eval_segments)
    
    # add no_eval segments to final_df
    no_eval = format_for_final_df(no_eval)
    no_eval_dedup = no_eval.drop_duplicates(subset=['RT_ID'], ignore_index=True)
    final_df = pd.concat([final_df, no_eval_dedup], ignore_index=True)
    
    #### WARNING: there is some overlap between final_df and the eval dfs, this will be taken care of after eval
    #### this is probably due to preserving some RT_IDs in final_df that overlap muliple genes some of which need
    #### to or can be evaluated and some that do not/cannot (i.e. genes too close together, etc.)
    
    # sort dataframes (required for count_bed)
    final_df = final_df.sort_values(by = ['RT_chr', 'RT_start']).rename(columns={'RT_chr': '#RT_chr'})
    eval_segments = eval_segments.sort_values(by = ['RT_chr', 'RT_start']).rename(columns={'RT_chr': '#RT_chr'})
    eval_start = eval_start.sort_values(by = ['chr', 'start']).rename(columns={'chr': '#chr'})
    eval_end = eval_end.sort_values(by = ['chr', 'start']).rename(columns={'chr': '#chr'})
    eval_gene = eval_gene.sort_values(by = ['chr', 'start']).rename(columns={'chr': '#chr'})
    
    # output files
    final_output = os.path.splitext(working_file)[0]+".noeval.bed"
    eval_output = os.path.splitext(working_file)[0]+".eval.bed"
    start_output = os.path.splitext(working_file)[0]+".starts.bed"
    end_output = os.path.splitext(working_file)[0]+".ends.bed"
    gene_output = os.path.splitext(working_file)[0]+".genes.bed"
    
    delete_if_exists(final_output)
    delete_if_exists(eval_output)
    delete_if_exists(start_output)
    delete_if_exists(end_output)
    delete_if_exists(gene_output)
    
    final_df.to_csv(final_output, sep='\t', header=True, index=False)
    eval_segments.to_csv(eval_output, sep='\t', header=True, index=False)
    eval_start.to_csv(start_output, sep='\t', header=True, index=False)
    eval_end.to_csv(end_output, sep='\t', header=True, index=False)
    eval_gene.to_csv(gene_output, sep='\t', header=True, index=False)


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('segment_file',
                        help='Path to RT_seg .bed file (expects .segmentRPKM.bed file).')
    parser.add_argument('-e', '--exon_expression_file',
                        help='Path to matrix containing 0h gene RPKM values, because for samples run with this modified script we do not have 6h data')
    parser.add_argument('sample',
                        help='String denoting the sample name (i.e. A6730h4001a).')
    
    args = parser.parse_args()

    main(args.segment_file, args.sample, args.exon_expression_file)
    
    # segments_file = None
    # exon_expression_file = None
    
    # main(segments_file, exon_expression_file)
    