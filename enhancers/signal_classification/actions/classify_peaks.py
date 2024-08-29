# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:52:13 2023

Description
-----------
This script defines the following functions which are used to classify eRNA uv_peaks:
    
    filter_list: filters list of files in dir for eRNA files with counts/rpkm
    filter_unipeak_list: filters list of files in dir for original eRNA unipeak files
    filter_allpeak_list: filters list of files in dir for original eRNA all peak summary files
    get_cl_list: creates a list of cell line names
    create_df_dicts: loads in files as pandas df and adds to dictionary
    stdev_and_mean: calculates mean and standard deviation of bipeak log10 rpkm ratio dist for all cell
                    lines and average values for all cell lines
    classify_bipeaks: classifies bipeaks by symmetry and standard deviation and adds columns to dfs
    classify_unipeaks: classifies unipeaks by rpkm ratio and bipeak StDev and adds columns to dfs
    fix_bipeaks: reformat df (column order) for final output
    fix_unipeaks: reformat df (merge back MACS2 columns removed and reorder columns) for final output
    fix_allpeaks: add classification columns to all peaks files
    update_files: runs functions to reformat dfs and creates new consensus matrices with class info
    save_from_dict: write updated cell line specific dfs to file from dictionaries
    
Peak Classifications
--------------------
    Peak type (assigned by previous scripts): single (AM), divergent, convergent, AcoversB, BcoversA (BM)
    Peak class: symmetrical, asymmetrical_plus, asymmetrical_minus, individual_plus, individual_minus
    Peak class confidence: unidir, bidir_ip/im_low3, bidir_ip/im_low, bidir_sy/ap/am, bidir_sy/ap/am_high1, bidir_sy/ap/am_high2, bidir_sy/ap/am_high3

@author: abmcs
"""

import argparse
import os
import pandas as pd
import re
from tqdm import tqdm
from statistics import mean

def filter_list(datalist):
    
    return [val for val in datalist if re.search(r'dELS_intersect.counts.rpkm.bed', val)]


def filter_unipeak_list(datalist):
    
    return [val for val in datalist if re.search(r'.unique_unipeaks_600.dELS_intersect.bed', val)]


def filter_allpeak_list(datalist):
    
    return [val for val in datalist if re.search(r'.unique_peaks_600.dELS_intersect.bed', val)]


def get_cl_list(bed_filepath_list):
    
    cellline_list = []
    for filepath in bed_filepath_list:
        cellline = filepath.rsplit(sep='\\', maxsplit=1)[1]
        cellline = cellline.split(sep='.')[0]
        cellline_list.append(cellline)
    cl_list = list(set(cellline_list))
    
    return cl_list


def create_df_dicts(fp_list, up_list, ap_list, celllines):
    
    bp_dict = {}
    up_dict = {}
    oup_dict = {}
    ap_dict = {}
    for cl in tqdm(celllines):
        # load in dfs for cell line
        for fp in fp_list:
            if cl not in fp:
                continue
            else:
                if 'bipeaks' in fp:
                    bp = pd.read_csv(fp, sep='\t', engine='python')
                    bp.rename(columns={'bipeak_class':'peak_type'}, inplace=True)
                if 'unipeaks' in fp:
                    up = pd.read_csv(fp, sep='\t', engine='python')
        bp_dict[cl] = bp
        up_dict[cl] = up
        # load in original unipeak files to create final df
        for up in up_list:
            if cl not in up:
                continue
            else:
                oup = pd.read_csv(up, sep='\t', engine='python')
        oup_dict[cl] = oup
        # load in all peak summary files to add classifications
        for p in ap_list:
            if cl not in p:
                continue
            else:
                ap = pd.read_csv(p, sep='\t', engine='python')
        ap_dict[cl] = ap
        
    return bp_dict, up_dict, oup_dict, ap_dict


def stdev_and_mean(df_dict):
    
    # get stdev for all dfs and add to list
    sd_list = []
    mean_list = []
    for df in df_dict.values():
        stdev = df['rpkm_ratio_log10'].std()
        sd_list.append(stdev)
        dist_mean = df['rpkm_ratio_log10'].mean()
        mean_list.append(dist_mean)
    # get mean across cell lines
    avg_sd = mean(sd_list)
    avg_mean = mean(mean_list)
    
    return sd_list, avg_sd, mean_list, avg_mean


def classify_bipeaks(bp_dict, stdev, dist_mean):
    
    # set positive cutoffs
    pos_cutoff = dist_mean + (stdev / 2)
    pos_high = dist_mean + (stdev / 4)
    pos_1sd = dist_mean + stdev
    pos_2sd = dist_mean + (stdev * 2)
    pos_3sd = dist_mean + (stdev * 3)
    # set negative cutoffs
    neg_cutoff = dist_mean - (stdev / 2)
    neg_high = dist_mean - (stdev / 4)
    neg_1sd = dist_mean - stdev
    neg_2sd = dist_mean - (stdev * 2)
    neg_3sd = dist_mean - (stdev * 3)
    # create updated df dict
    updated_dict = {}
    for cl, df_original in bp_dict.items():
        df = df_original.copy()
        df['peak_class'] = ''
        df['peak_class_confidence'] = ''
        for idx, row_val in tqdm(enumerate(df.iterrows()), total=len(df), desc='Classifying peaks:'+cl):
            row = row_val[1]
            if row['rpkm_ratio_log10'] > pos_cutoff:
                df.loc[idx, 'peak_class'] = 'asymmetrical_plus'
                if row['rpkm_ratio_log10'] < pos_1sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_ap'
                if row['rpkm_ratio_log10'] > pos_1sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_ap_high1'
                if row['rpkm_ratio_log10'] > pos_2sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_ap_high2'
                if row['rpkm_ratio_log10'] > pos_3sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_ap_high3'
            if row['rpkm_ratio_log10'] < neg_cutoff:
                df.loc[idx, 'peak_class'] = 'asymmetrical_minus'
                if row['rpkm_ratio_log10'] > neg_1sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_am'
                if row['rpkm_ratio_log10'] < neg_1sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_am_high1'
                if row['rpkm_ratio_log10'] < neg_2sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_am_high2'
                if row['rpkm_ratio_log10'] < neg_3sd:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_am_high3'
            if row['rpkm_ratio_log10'] < pos_cutoff and row['rpkm_ratio_log10'] > neg_cutoff:
                df.loc[idx, 'peak_class'] = 'symmetrical'
                if row['rpkm_ratio_log10'] < pos_high and row['rpkm_ratio_log10'] > neg_high:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_sy_high'
                else:
                    df.loc[idx, 'peak_class_confidence'] = 'bidir_sy'
        updated_dict[cl] = df
    
    return updated_dict


def classify_unipeaks(up_dict, bidir_stdev, bidir_dist_mean):
    # get bidirectional cutoff (simply for classification)
    neg_3sd = bidir_dist_mean - (bidir_stdev * 3)
    # create updated df dict
    updated_dict = {}
    for cl, df_original in up_dict.items():
        df = df_original.copy()
        df['peak_class'] = ''
        df['peak_class_confidence'] = ''
        for idx, row_val in tqdm(enumerate(df.iterrows()), total=len(df), desc='Classifying unipeaks:'+cl):
            row = row_val[1]
            if row['orig_unipeak_strand'] == '+':
                df.loc[idx, 'peak_class'] = 'individual_plus'
                if row['rpkm_ratio'] == 0:
                    df.loc[idx, 'peak_class_confidence'] = 'unidir_ip'
                else:
                    if row['rpkm_ratio_log10'] < neg_3sd:
                        df.loc[idx, 'peak_class_confidence'] = 'bidir_ip_low3'
                    if row['rpkm_ratio_log10'] > neg_3sd:
                        df.loc[idx, 'peak_class_confidence'] = 'bidir_ip_low'
            if row['orig_unipeak_strand'] == '-':
                df.loc[idx, 'peak_class'] = 'individual_minus'
                if row['rpkm_ratio'] == 0:
                    df.loc[idx, 'peak_class_confidence'] = 'unidir_im'
                else:
                    if row['rpkm_ratio_log10'] < neg_3sd:
                        df.loc[idx, 'peak_class_confidence'] = 'bidir_im_low3'
                    if row['rpkm_ratio_log10'] > neg_3sd:
                        df.loc[idx, 'peak_class_confidence'] = 'bidir_im_low'
        updated_dict[cl] = df
    
    return updated_dict


def fix_bipeaks(bp_dict, cl):
    
    bp_df = bp_dict[cl]
    bp_df = bp_df.sort_values(by=['bipeak_chr','bipeak_start'], ignore_index=True)
    bp_df_rn = bp_df.rename(columns={'counts_peak1_counts':'peak1_counts', 'counts_peak1_density':'peak1_density',
                                     'counts_peak1_RPKM':'peak1_RPKM', 'counts_peak2_counts':'peak2_counts',
                                     'counts_peak2_density':'peak2_density', 'counts_peak2_RPKM':'peak2_RPKM'})
    final_bp = bp_df_rn[['bipeak_chr', 'bipeak_start', 'bipeak_end', 'bipeak_id', 'bipeak_score', 'bipeak_strand', 'bipeak_midpoint',
           'peak1_chr', 'peak1_start', 'peak1_end', 'peak1_id', 'peak1_score', 'peak1_strand', 'peak1_5p_start', 'peak1_5p_end',
           'peak1_signal_value', 'peak1_pval', 'peak1_qval', 'peak1_point_source', 'peak1_counts', 'peak1_density', 'peak1_RPKM',
           'peak2_chr', 'peak2_start', 'peak2_end', 'peak2_id', 'peak2_score', 'peak2_strand', 'peak2_5p_start', 'peak2_5p_end',
           'peak2_signal_value', 'peak2_pval', 'peak2_qval', 'peak2_point_source', 'peak2_counts', 'peak2_density', 'peak2_RPKM',
           'distance', 'distance_signed', 'peak_type', 'peak_id', 'rpkm_ratio', 'rpkm_ratio_log10', 'peak_class', 'peak_class_confidence']]
    
    return final_bp


def fix_unipeaks(up_dict, oup_dict, cl):
    
    up_df = up_dict[cl]
    up_df = up_df.sort_values(by=['div_unipeak_chr','div_unipeak_start'], ignore_index=True)
    oup_df = oup_dict[cl]
    # fix column names for merge
    up_df.columns = up_df.columns.str.replace('orig_unipeak_|orig_peak_|orig_|unipeaks_peak_', 'peak_', regex=True)
    up_df.columns = up_df.columns.str.replace('div_unipeak_', 'div_', regex=True)
    up_df_rn = up_df.rename(columns={'unipeak_class':'peak_type'})
    oup_df_rn = oup_df.rename(columns={'X5p_start':'peak_5p_start', 'X5p_end':'peak_5p_end', 'unipeak_class':'peak_type'})
    # merge dfs
    merged_up = oup_df_rn.merge(up_df_rn, how='inner', on=['peak_id', 'peak_5p_start', 'peak_5p_end'])
    # drop redundant columns, update column values, and reorder columns
    subset_up = merged_up.drop(columns=['peak_start', 'peak_end', 'peak_strand'])
    subset_up['div_score'] = 0
    final_up = subset_up[['unipeak_chr', 'unipeak_start', 'unipeak_end', 'unipeak_id', 'unipeak_score', 'unipeak_strand',
                          'peak_id', 'peak_5p_start', 'peak_5p_end', 'signal_value', 'pval', 'qval', 'point_source',
                          'peak_type', 'peak_counts', 'peak_density', 'peak_RPKM',
                          'div_chr', 'div_start', 'div_end', 'div_peak_id', 'div_score', 'div_strand', 'div_5p_start',
                          'div_5p_end', 'div_counts', 'div_density', 'div_RPKM', 
                          'rpkm_ratio', 'rpkm_ratio_log10', 'peak_class', 'peak_class_confidence']]
    
    return final_up


def fix_allpeaks(ap_dict, fixed_up_dict, fixed_bp_dict, cl):
    
    ap_df = ap_dict[cl]
    up_df = fixed_up_dict[cl].copy()
    bp_df = fixed_bp_dict[cl].copy()
    # subset unipeak and bipeak dfs to get required columns
    up_df_sub = up_df[['unipeak_chr', 'unipeak_start', 'unipeak_end', 'unipeak_id', 'unipeak_score',
                       'unipeak_strand', 'peak_id', 'peak_type', 'peak_class', 'peak_class_confidence']].copy()
    bp_df_sub = bp_df[['bipeak_chr', 'bipeak_start', 'bipeak_end', 'bipeak_id', 'bipeak_score',
                       'bipeak_strand', 'peak_id', 'peak_type', 'peak_class', 'peak_class_confidence']].copy()
    # rename columns
    ap_df_rn = ap_df.rename(columns={'class':'peak_type'})
    up_df_rn = up_df_sub.rename(columns={'unipeak_id':'class_id'})
    up_df_rn.columns = up_df_rn.columns.str.replace('unipeak_', '', regex=True)
    bp_df_rn = bp_df_sub.rename(columns={'bipeak_id':'class_id'})
    bp_df_rn.columns = bp_df_rn.columns.str.replace('bipeak_', '', regex=True)
    bp_df_rn.rename(columns={'midpoint':'bipeak_midpoint'}, inplace=True)
    # merge unipeak and bipeak dfs with all peak to add class columns
    cat_df = pd.concat([bp_df_rn, up_df_rn], ignore_index=True)
    merged_ap = ap_df_rn.merge(cat_df, how='left', on=['chr', 'start', 'end', 'class_id', 'score', 'strand', 'peak_id', 'peak_type'])
    
    return merged_ap


def update_files(bp_dict, up_dict, oup_dict, ap_dict, matrix, celllines):
    
    # fix bipeak files (column order)
    fixed_bp_dict = {}
    for cl in tqdm(celllines, desc='Fixing bipeak dfs'):
        final_bp = fix_bipeaks(bp_dict, cl)
        fixed_bp_dict[cl] = final_bp
    
    # fix unipeak files with divergent ratios (start with original unipeak files and add stuff from Karan's files)
    fixed_up_dict = {}
    for cl in tqdm(celllines, desc='Fixing unipeak dfs'):
        final_up = fix_unipeaks(up_dict, oup_dict, cl)
        fixed_up_dict[cl] = final_up
    
    # fix all peak files (add classifications)
    fixed_ap_dict = {}
    for cl in tqdm(celllines, desc='Fixing all peak dfs'):
        final_ap = fix_allpeaks(ap_dict, fixed_up_dict, fixed_bp_dict, cl)
        fixed_ap_dict[cl] = final_ap
    
    # add class and class confidence to consensus matrix
    type_matrix = matrix.copy()
    class_matrix = matrix.copy()
    conf_matrix = matrix.copy()
    for cl in tqdm(celllines, desc='Updating consensus matrix'):
        ap_df = fixed_ap_dict[cl].copy()
        for idx, row_val in tqdm(enumerate(ap_df.iterrows()), total=len(ap_df), desc=cl):
            row = row_val[1]
            peak = row['class_id']
            mat_loc = matrix.loc[matrix[cl] == peak, cl]
            if len(mat_loc) == 0:
               type_matrix.replace(fr'(?<![^,]){peak}(?![^,])', row['peak_type'], regex=True, inplace=True)
               class_matrix.replace(fr'(?<![^,]){peak}(?![^,])', row['peak_class'], regex=True, inplace=True)
               conf_matrix.replace(fr'(?<![^,]){peak}(?![^,])', row['peak_class_confidence'], regex=True, inplace=True)
            else:
                type_matrix.loc[type_matrix[cl] == peak, cl] = row['peak_type']
                class_matrix.loc[class_matrix[cl] == peak, cl] = row['peak_class']
                conf_matrix.loc[conf_matrix[cl] == peak, cl] = row['peak_class_confidence']
                
    return fixed_bp_dict, fixed_up_dict, fixed_ap_dict, type_matrix, class_matrix, conf_matrix


def save_from_dict(df_dict, suffix):
    
    for cl, df in df_dict.items():
        filename = 'C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/' + cl + suffix
        df.to_csv(filename, sep='\t', header=True, index=False)
    
    return


def main(fp_list, up_list, ap_list):
    
    # load consensus region matrix
    matrix = pd.read_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/peak_overlaps/BruUV_eRNA_16_cell_lines.bed', sep='\t', engine='python')
    
    celllines = get_cl_list(fp_list)
    
    bipeaks, unipeaks, original_unipeaks, all_peaks = create_df_dicts(fp_list, up_list, ap_list, celllines)
    stdev_list, stdev, mean_list, dist_mean = stdev_and_mean(bipeaks)
    
    # make bipeak classifications based on half of the mean stdev (stdev) from the average mean of all distributions (mean)
    bipeaks_class = classify_bipeaks(bipeaks, stdev, dist_mean)
    
    # make unipeak classifications based on rpkm ratio of 0 and then call rest low confidence bipeaks
    unipeaks_class = classify_unipeaks(unipeaks, stdev, dist_mean)
    
    # fix/update files with classifications
    final_bipeaks, final_unipeaks, final_allpeaks, type_matrix, class_matrix, conf_matrix = update_files(bipeaks_class, unipeaks_class, original_unipeaks, all_peaks, matrix, celllines)
    
    # output all updated files with classifications
    save_from_dict(final_bipeaks, '.unique_bipeaks_600.dELS_intersect.counts.rpkm.class.bed')
    save_from_dict(final_unipeaks, '.unique_unipeaks_600.dELS_intersect.counts.rpkm.class.bed')
    save_from_dict(final_allpeaks, '.unique_peaks_600.dELS_intersect.counts.rpkm.class.bed')
    
    # output matrices
    type_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakType.bed', sep='\t', header=True, index=False, na_rep='NA')
    type_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakType.tsv', sep='\t', header=True, index=False, na_rep='NA')
    class_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakClass.bed', sep='\t', header=True, index=False, na_rep='NA')
    class_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakClass.tsv', sep='\t', header=True, index=False, na_rep='NA')
    conf_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakClassConfidence.bed', sep='\t', header=True, index=False, na_rep='NA')
    conf_matrix.to_csv('C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/final_files/BruUV_eRNA_16_cell_lines_peakClassConfidence.tsv', sep='\t', header=True, index=False, na_rep='NA')
    
    return


if __name__=="__main__":
    
    # parser = argparse.ArgumentParser()
    # parser.add_argument('directory',
                        # help='Path to original file directory.')
    
    args = parser.parse_args()
    
    file_dir = "C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/ratios"
    fp_list_all = [os.path.abspath(os.path.join(file_dir, p)) for p in os.listdir(file_dir)]
    fp_list = filter_list(fp_list_all)
    
    ip_dir = "C:/Users/abmcs/Desktop/shared_biosamples/UV stuff/bipeak/intersected_peak_files"
    up_list_all = [os.path.abspath(os.path.join(ip_dir, p)) for p in os.listdir(ip_dir)]
    up_list = filter_unipeak_list(up_list_all)
    
    ap_list = filter_allpeak_list(up_list_all)
    
    
    main(fp_list, up_list, ap_list)