# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:14:08 2022

Description
-----------
This script defines the following functions which are used to... :
    
    load_seg_files: loads and formats segment files needed for this script
    load_exp_files: loads and formats expression files needed for this script
    load_close_files: load and subsets files with RT distance to closest downstream gene
    load_rev_files: loads and subsets files with RT intersections with genes on the opposite strand
    delete_if_exists: deletes file if already exists in working directory
    add_exp_col: adds column containing 6h expression data
    update_classes: reassigns RT_classes for all rows
    merge_duplicate_segs: propagates RT_classes to any matching RT_segs and deduplicates df
    
@author: abmcs
"""


import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
from collections import Counter
import warnings


def load_seg_files(filepath):
    
    '''
    loads in pandas dataframes of each of the segment files required for this step, adds a header, adds 2 columns
    containing 1.) the gene_length (bp) of the overlapping genes found, and 2.) containing the number of gene 
    overlaps for each RT_segment, and restructures the column orgnaization of the df.
    '''
    
    header = ['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length', 'old_counts', 'old_density', 'old_rpkm', 'original_RTseg',
              'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
              'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_length', 'X',
              'seglen', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm']
    
    df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    
    # create column with overlapping gene length
    df['ol_gene_length'] = pd.to_numeric(df['ol_end']) - pd.to_numeric(df['ol_start'])
    
    # create column with counts of gene overlaps per segment (segs with no overlap will say 1, fixed later)
    counts = df['RT_ID'].value_counts().to_frame()
    counts.reset_index(inplace=True)
    counts.rename(columns={'RT_ID':'num_geneOL', 'index':'RT_ID'}, inplace=True)
    df = pd.merge(df, counts, on=['RT_ID'], how='outer')
    
    # set RT_length as value from count_bed in case I messed up somewhere
    df['RT_length'] = df['seglen']
    
    formatted_df = df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm', 'original_RTseg',
                       'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
                       'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_gene_length', 'ol_length', 'num_geneOL']]
    
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


def load_close_files(filepath):
    
    '''
    loads in pandas dataframes of the files containing the closest downstream genes to the fixed RT_segments, and
    subsets dataframe to contain only relevant columns for look-up (used to refine class I annotation); all zero
    distances are overlaps
    '''
    
    header = ['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length', 'old_counts', 'old_density', 'old_rpkm', 'original_RTseg',
              'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
              'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_length', 'X',
              'seglen', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm',
              'ds_chr', 'ds_start', 'ds_end', 'ds_gene', 'ds_score', 'ds_strand', 'ds_feat_type', 'ds_dist']
    
    df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    
    formatted_df = df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length',
                       'ds_chr', 'ds_start', 'ds_end', 'ds_gene', 'ds_score', 'ds_strand', 'ds_feat_type', 'ds_dist']]
    
    output_df = formatted_df.drop_duplicates(subset=['RT_ID'], ignore_index=True)
    
    return output_df


def load_rev_files(filepath):
    
    '''
    loads in pandas dataframes of the files containing genes on the opposite strand that overlap the fixed 
    RT_segments, and subsets dataframe to contain only relevant columns for classification (class IV)
    '''
    
    header = ['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length', 'old_counts', 'old_density', 'old_rpkm', 'original_RTseg',
              'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm',
              'ol_chr', 'ol_start', 'ol_end', 'ol_gene_ID', 'ol_score', 'ol_strand', 'ol_feat_type', 'ol_length', 'X',
              'seglen', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm',
              'rv_chr', 'rv_start', 'rv_end', 'rv_gene', 'rv_score', 'rv_strand', 'rv_feat_type', 'rv_len']
    
    df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
    
    formatted_df = df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_class', 'RT_length',
                       'rv_chr', 'rv_start', 'rv_end', 'rv_gene', 'rv_score', 'rv_strand', 'rv_feat_type', 'rv_len']]
    
    return formatted_df


def delete_if_exists(filepath):
    
    '''
    removes existing output file if filename matches one in directory
    '''
    
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    
    return False


def add_exp_col(df, df_exon_exp):
    
    '''
    Parameters
    ----------
    df: pandas df with updated segments, segment counts/rpkm, and gene overlap info
    df_exon_exp: pandas df containing BruChase exon RPKM values for evaluating gene expression at 6h timepoint
    
    Description
    -----------
    This function checks genes overlapping new RT_segments for expression in 6h data (rpkm > 0.5) and adds binary
    column describing expression (yes | no).
            
    Variables
    ---------
    6h expression cutoff can be changed (currently >0.5) --> modified to 0.25 becuase this script uses 0h gene rpkm
    '''
    
    df = df.reset_index(drop=True)
    
    # add expressed column if needed
    if 'expressed' not in df:
        df['expressed'] = np.nan
    
    expressed_df_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Check for expression in 6h data"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        if row['ol_chr'] == '.':
            expressed_df_rows.append(row)
            continue
        gene_name = df.loc[idx, 'ol_gene_ID'].split('/')[0]
        gene_id = df.loc[idx, 'ol_gene_ID'].split('/')[1]
        exp = df_exon_exp.loc[((df_exon_exp['gene_name'] == gene_name) & (df_exon_exp['gene_ID'] == gene_id)), 'rpkm']
        # make sure there is a gene exp value, if not label as missing
        if exp.empty:
            row['expressed'] = 'missing'
            expressed_df_rows.append(row)
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
    
    return expression_df


def add_rev_exp_col(df_rev, df_gene_exp):
    
    '''
    Parameters
    ----------
    df_rev: pandas df with information about genes overlapping RT_segments on the opposite strand
    df_gene_exp: pandas df containing RPKM values for evaluating gene expression at 0h timepoint
    
    Description
    -----------
    This function checks genes overlapping new RT_segments on the opposite strand for expression in 0h data
    (rpkm > 0.25) and adds binary column describing expression (yes | no).
            
    Variables
    ---------
    0h expression cutoff can be changed (currently >0.25)
    '''
    
    df = df_rev.copy()
    
    # add expressed column if needed
    if 'expressed' not in df:
        df['expressed'] = np.nan
    
    expressed_df_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Check for expression in 0h data"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        if row['rv_chr'] == '.':
            expressed_df_rows.append(row)
            continue
        gene_name = df.loc[idx, 'rv_gene'].split('/')[0]
        gene_id = df.loc[idx, 'rv_gene'].split('/')[1]
        exp = df_gene_exp.loc[((df_gene_exp['gene_name'] == gene_name) & (df_gene_exp['gene_ID'] == gene_id)), 'rpkm']
        # make sure there is a gene exp value, if not label as missing
        if exp.empty:
            row['expressed'] = 'missing'
            expressed_df_rows.append(row)
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
    
    return expression_df


def reverse_classification(df):
    
    #df = df_rev_exp.copy()
    df['rv_class'] = np.nan
    class_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Classify segments by reverse gene"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        if row['rv_chr'] == '.':
            continue
        # IVa --> rv_gene is short
        if row['rv_len'] < 1000:
            row['rv_class'] = 'IVa'
            class_rows.append(row)
            continue
        # IVb --> rv_gene not expressed and any overlap
        if row['expressed'] == 'no':
            row['rv_class'] = 'IVb'
        
        if row['expressed'] == 'yes': 
            if row['rv_strand'] == '+':
                # IVc --> rv_gene expressed, TSS does not overlap RT_seg, and is >2kb from RT_seg end
                if row['rv_start'] < int(row['RT_start']) - 2000:
                    row['rv_class'] = 'IVc'
                # IVd --> rv_gene expressed, TSS does not overlap RT_seg, and is <2kb from RT_seg end
                if (row['rv_start'] < row['RT_start']) and not (row['rv_start'] < int(row['RT_start']) - 2000):
                    row['rv_class'] = 'IVd'
                # IVe --> rv_gene expressed and TSS overlaps RT_seg (full overlap)
                if (row['rv_start'] >= row['RT_start']) and (row['rv_end'] <= row['RT_end']):
                    row['rv_class'] = 'IVe'
                # IVf --> rv_gene expressed and TSS overlaps RT_seg (partial overlap)
                if (row['rv_start'] >= row['RT_start']) and (row['rv_end'] > row['RT_end']):
                    row['rv_class'] = 'IVf'
            else:
                # IVc --> rv_gene expressed, TSS does not overlap RT_seg, and is >2kb from RT_seg end
                if row['rv_end'] > int(row['RT_end']) + 2000:
                    row['rv_class'] = 'IVc'
                # IVd --> rv_gene expressed, TSS does not overlap RT_seg, and is <2kb from RT_seg end
                if (row['rv_end'] > row['RT_end']) and not (row['rv_end'] > int(row['RT_end']) + 2000):
                    row['rv_class'] = 'IVd'
                # IVe --> rv_gene expressed and TSS overlaps RT_seg (full overlap)
                if (row['rv_end'] <= row['RT_end']) and (row['rv_start'] >= row['RT_start']):
                    row['rv_class'] = 'IVe'
                # IVf --> rv_gene expressed and TSS overlaps RT_seg (partial overlap)
                if (row['rv_end'] <= row['RT_end']) and (row['rv_start'] < row['RT_start']):
                    row['rv_class'] = 'IVf'
        # check if expression is missing and add as class V  
        if row['expressed'] == 'missing': 
            row['rv_class'] = 'V'
        class_rows.append(row)
        
    class_df = pd.DataFrame(class_rows)
    
    # check that all reverse classes in all rows are the same
    id_set = set(class_df['RT_ID'])
    for seg in id_set:
        # seg = 'chr1:766250-826205'
        df_sub = class_df[class_df['RT_ID'] == seg]
        df_sub = df_sub.reset_index(drop=True)
        if all(df_sub['rv_class'] == df_sub['rv_class'][0]):
            continue
        else:
            classes = ','.join(sorted(list(df_sub['rv_class'].unique())))
            class_df.loc[class_df['RT_ID'] == seg, 'rv_class'] = classes
    
    # make lookup dict (key = RT_ID, value = classes)
    # class_df = df_rev_class.copy()
    class_df_sub = class_df[['RT_ID', 'rv_class']].drop_duplicates(ignore_index=True)
    class_dict = dict(zip(class_df_sub.RT_ID, class_df_sub.rv_class))
    
    return class_dict


def update_classes(df, df_close, rev_dict):
    
    '''
    Parameters
    ----------
    df: pandas df containing updated segment coordinates and overlapping gene information
    df_close: df containing distance to closest downstream gene (for class I designations)
    rev_dict: dictionary where keys are RT_segs and values are rev strand classifications (from reverse_classification())
    
    Description
    -----------
    This function adds a new_RT_class column and reclassifies all segments according to established rules (see 
    comments below and/or Read_through Classes google slides image). Redundant classifications are then removed.
        
    Variables
    ---------
    Gene length and overlap cutoffs for classifications can be changed here (small genes <1kb; small overlaps <2kb)
    '''
    
    df['new_RT_class'] = np.nan
    # reclassify segments
    classified_rows = []
    # for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Re-classify segments"):
    for row_ind, row in df.iterrows():
        # row = row_vals[1]
        # fix number of gene overlaps for segs that do not overlap a gene and classify according to closest ds gene
        if row['ol_chr'] == '.':
            row['num_geneOL'] = 0
            if df_close.loc[df_close['RT_ID'] == row['RT_ID'], 'ds_dist'].item() > 5000:
                if pd.isna(row['new_RT_class']):
                    row['new_RT_class'] = 'Ia'
                else:
                    row['new_RT_class'] = row['new_RT_class'] + ',Ia'
            else:
                if pd.isna(row['new_RT_class']):
                    row['new_RT_class'] = 'Ib'
                else:
                    row['new_RT_class'] = row['new_RT_class'] + ',Ib'
        
        # add segment classifications from reverse strand overlap
        if row['RT_ID'] in rev_dict:
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = rev_dict[row['RT_ID']]
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',' + rev_dict[row['RT_ID']]
        
        # classify segments that were re-evaluated
        if any(x in row['RT_class'] for x in ['IId','X']):
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IId'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IId'
        
        # classify small genes and do not further classify these rows
        if row['ol_gene_length'] < 1000 and not row['ol_chr'] == '.':
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIa'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIa'
            classified_rows.append(row)
            continue
        # classify partial overlaps of not expressed genes
        if (row['ol_gene_length'] != row['ol_length']) and (row['expressed'] == 'no'):
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIb'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIb'
        # classify full overlaps with not expressed genes
        if (row['ol_gene_length'] == row['ol_length']) and (row['expressed'] == 'no'):
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIc'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIc'
        
        ### IId moved above to ensure all reeval segs are captured
        
        # classify segments that have a small overlap with an expressed gene
        if (row['ol_length'] < 2000) and (row['expressed'] == 'yes') and (row['ol_gene_length'] > 2000):
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIe'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIe'
        # # classify partial overlaps of genes that do not have expression at the start of the gene
        # if 'IIf' in row['RT_class']:
        #     if pd.isna(row['new_RT_class']):
        #         row['new_RT_class'] = 'IIf'
        #     else:
        #         row['new_RT_class'] = row['new_RT_class'] + ',IIf'
        # classify ambiguous segment where overlapping gene is expressed
        if (row['expressed'] == 'yes') and not any(x in row['RT_class'] for x in ['IId','X']):
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIIa'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIIa'
        # classify ambiguous segments with overlapping genes
        if (row['expressed'] == 'yes') and not any(x in row['RT_class'] for x in ['IId','X']):
            # check if gene is completely overlapped by ol_gene
            if row['gene_start'] >= row['ol_start'] and (row['gene_end'] <= row['ol_end']):
                if pd.isna(row['new_RT_class']):
                    row['new_RT_class'] = 'IIIb'
                else:
                    row['new_RT_class'] = row['new_RT_class'] + ',IIIb'
            if row['RT_strand'] == '+':
                # check if gene is partially overlapped by ol_gene
                if (row['gene_start'] < row['ol_start']) and (row['gene_end'] > row['ol_start']):
                    if pd.isna(row['new_RT_class']):
                        row['new_RT_class'] = 'IIIc'
                    else:
                        row['new_RT_class'] = row['new_RT_class'] + ',IIIc'
            else:
                # check if gene is partially overlapped by ol_gene
                if (row['gene_end'] > row['ol_end']) and (row['gene_start'] < row['ol_end']):
                    if pd.isna(row['new_RT_class']):
                        row['new_RT_class'] = 'IIIc'
                    else:
                        row['new_RT_class'] = row['new_RT_class'] + ',IIIc'
        # check to make sure overlapping genes don't match the RT gene, if yes classify as 'IIId'
        if row['gene_name'] == row['ol_gene_ID'].split('/')[0]:
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'IIId'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',IIId'
        # assign segments overlapping genes with no expression info to class 'IIIe' (changed 09072023, class V now describes missing expression info; should not exist anymore)
        if row['expressed'] == 'missing':
            if pd.isna(row['new_RT_class']):
                row['new_RT_class'] = 'V'
            else:
                row['new_RT_class'] = row['new_RT_class'] + ',V'
        classified_rows.append(row)
            
    classified_df = pd.DataFrame(classified_rows)
    
    # remove redundancies
    concise_row = []
    # remove redundant classes
    # for row_ind, row_vals in tqdm(enumerate(classified_df.iterrows()), total=len(classified_df), desc="Remove redundancies"):
    for row_ind, row in classified_df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        # check for duplicates and remove them
        string = classified_df.loc[idx, 'new_RT_class']
        new_classes = sorted(list(set(string.split(','))))
        # if contains IId and IIb, remove IIb (should not occur)
        if all(x in new_classes for x in ['IIb', 'IId']):
                new_classes.remove('IIb')
        # if contains IIb and IIe/IIf, only needs IIe/IIf (should not occur)
        if 'IIb' in new_classes and any(x in new_classes for x in ['IIe', 'IIf']):
            new_classes.remove('IIb')
        # if it is IIe/IIf/IIIb/IIIc, no need for IIIa
        if 'IIIa' in new_classes and any(x in new_classes for x in ['IIe', 'IIf', 'IIIb', 'IIIc']):
            new_classes.remove('IIIa')
        # if it is IIIc, it's not IIe (this means overlap is at end of gene not start)
        if 'IIe' in new_classes and any(x in new_classes for x in ['IIIb', 'IIIc']):
        #if any((x == 'IIIc' and x == 'IIe') for x in new_classes):
            new_classes.remove('IIe')
        # change all class IIIc/d to class IIIb (new as of 09072023, generalizing overlapping gene classes)
        if 'IIIc' in new_classes:
            new_classes = list(map(lambda x: x.replace('IIIc', 'IIIb'), new_classes))
        if 'IIId' in new_classes:
            new_classes = list(map(lambda x: x.replace('IIId', 'IIIb'), new_classes))
        final_classes = list(set(new_classes))
        # join back together and set as RT_class
        new_string = ','.join(sorted(final_classes, key=str.casefold))
        row['new_RT_class'] = new_string
        concise_row.append(row)
    concise_df = pd.DataFrame(concise_row)
    
    updated_df = concise_df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'new_RT_class', 'RT_length', 'RTseg_counts', 'RTseg_density', 'RTseg_rpkm', 'original_RTseg',
                             'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm', 'num_geneOL']]
    
    final_df = updated_df.rename(columns={"new_RT_class":"RT_class"})
    
    return final_df


def merge_duplicate_segs(df):
    
    '''
    Parameters
    ----------
    df: segments df with updated RT_class column (output from update_classes())
    
    Description
    -----------
    This function checks segments with duplicate entries (overlap multiple genes) and assigns all classes to all
    rows, then removes duplicate rows and outputs this as final df.
        
    Variables
    ---------
    No variables to change in this function.
    '''
    
    # make sure that all new RT_IDs have the same classes in all rows
    new_ids = set(df['RT_ID'])
    for seg in new_ids:
        df_sub = df[df['RT_ID'] == seg]
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
            df.loc[df['RT_ID'] == seg, 'RT_class'] = RT_class
    
    final_df = df.drop_duplicates(ignore_index=True)
    
    return final_df


def main(segments_file, sample, gene_exp_file, closest_gene_file, reverse_gene_file, exon_exp_file=None):
    
    '''
    Parameters
    ----------
    filepath: [arg1] path to RT_seg .bed file with gene overlap and count_bed info
              (expects .segmentRPKM.fixedRPKM.bed file)
    closest_gene_file: [arg2]
    reverse_gene_file: [arg3]
    exon_exp_file: [arg4] path to file containing 6h exon RPKM values
              (expects .featureTypes.bed file)
    
    Description
    -----------
    This function updates RT_segments, reapplies cutoffs and reclassifies all segments. Classifications are
    checked for redundancy and propagated to all matching rows. Dataframe is then deduplicated and output as final
    iteration of read-through segments.
    
    Variables
    ---------
    
    '''
    
    if segments_file is None:
        segments_file = "~/Desktop/shared_biosamples/read-through/fused_dec/HCT1160h4001a.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed"
    if gene_exp_file is None:
        gene_exp_file = "~/Desktop/shared_biosamples/ALLrpkm_ENCODE16CL_0h_genes.bed"
    if closest_gene_file is None:
        closest_gene_file = "~/Desktop/shared_biosamples/read-through/fused_dec/HCT1160h4001a.RTsegments.TESintersect.segmentRPKM.closestDS.bed"
    if reverse_gene_file is None:
        reverse_gene_file = "~/Desktop/shared_biosamples/read-through/fused_dec/HCT1160h4001a.RTsegments.TESintersect.segmentRPKM.reverseOL.bed"
    # build in workaround for no 6h data
    if exon_exp_file is None:
        exon_exp_file = gene_exp_file
        assay = None
    else:
        assay = '6h'
    
    
    
    df_seg = load_seg_files(segments_file)
    df_exon_exp = load_exp_files(exon_exp_file, sample, assay=assay)
    df_gene_exp = load_exp_files(gene_exp_file, sample)
    df_close = load_close_files(closest_gene_file)
    df_rev = load_rev_files(reverse_gene_file)
    
    # add gene exp info to reverse strand df and assign reverse classes
    df_rev_exp = add_rev_exp_col(df_rev, df_gene_exp)
    rev_dict = reverse_classification(df_rev_exp)
    
    # take RTseg count and gene rpkm cutoff again, in case any new segments don't make the cut
    df_cutoff = df_seg[(df_seg['RTseg_counts'] > 10) & (df_seg['gene_rpkm'] > 0.25)]
    
    # add gene exp info for overlapping genes
    df_exp = add_exp_col(df_cutoff, df_exon_exp)
    
    # do some checks to make sure all info lines up and redefine RT_classes
    updated_df = update_classes(df_exp, df_close, rev_dict)
    
    # remove duplicate segments
    final_df = merge_duplicate_segs(updated_df)
    
    # sort dataframe
    final_df = final_df.sort_values(by = ['RT_chr', 'RT_start']).rename(columns={"RT_chr":"#RT_chr"})
    
    
    # TEST: check that number of classifications does not exceed the number of overlapping genes (except class I)
    # will not be true for some reeval segs, so will need to reclassify them
    # IId, IIIa, IIIb, IIIc, IVa-e are all extra classifications
    error_row = []
    # for row_ind, row_vals in tqdm(enumerate(final_df.iterrows()), total=len(final_df), desc="Final check"):
    for row_ind, row in final_df.iterrows():
        # row = row_vals[1]
        idx = row_ind
        string = final_df.loc[idx, 'RT_class']
        new_classes = sorted(list(set(string.split(','))))
        test = Counter(new_classes)
        classes = ['I', 'IIa', 'IIb', 'IIc', 'IIe', 'IIf']
        count = sum(test[x] for x in classes)
        if row['num_geneOL'] < count and not 'I' in new_classes:
            error_row.append(row)
            warnings.warn('The number of RT_classes exceeds the number of overlapping genes. See error_rows in final_filter_and_classification.py.')
    
    # output file
    output_string = os.path.splitext(segments_file)[0]
    final_output = output_string.rsplit(sep='.', maxsplit=1)[0]+".final.bed"
     
    delete_if_exists(final_output)
    
    final_df.to_csv(final_output, sep='\t', header=True, index=False)
    
    
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('segment_file',
                        help='Path to RT_seg .bed file (expects .segmentRPKM.fixedRPKM.bed file).')
    parser.add_argument('sample',
                        help='String denoting the sample name (i.e. A6730h4001a).')
    parser.add_argument('-e', '--exon_expression_file',
                        help='Path to matrix containing 6h exon RPKM values.')
    parser.add_argument('gene_expression_file',
                        help='Path to matrix containing 0h gene RPKM values.')
    parser.add_argument('closest_gene_file',
                        help='Path to file containing closest genes downstream of RT_seg (expects  file).')
    parser.add_argument('reverse_gene_file',
                        help='Path to file containing genes overlapping RT_segs on the opposite strand (expects  file).')

    
    args = parser.parse_args()

    main(args.segment_file, args.sample, args.gene_expression_file, args.closest_gene_file, args.reverse_gene_file, args.exon_expression_file)
    
    # segments_file = None
    # exon_expression_file = None
    # gene_expression_file = None
    # closest_gene_file = None
    # reverse_gene_file = None
    
    # main(segments_file, exon_expression_file, gene_expression_file, closest_gene_file, reverse_gene_file)
                