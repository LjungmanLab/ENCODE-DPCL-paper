# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 09:30:46 2022

Filter bi-peaks

@author: abmcs
"""

# load data as pandas dataframe

import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np

def delete_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    return False


def load_from_bed(bed_filepath, header):
    
    df = pd.read_csv(bed_filepath, sep='\t', names=header, engine='python')
    rep_dict = dict.fromkeys(header[:5] + header[6:], '.')
    df = df.replace(rep_dict, value=np.NaN)
    numeric_cols = ['bipeak_start', 'bipeak_end', 'bipeak_midpoint', 'peak1_5p_start', 'peak1_5p_end', 'peak1_score', 'peak1_signal_value', 'peak1_pval', 'peak1_qval', 'peak1_point_source', 'peak1_start', 'peak1_end', 'peak2_5p_start', 'peak2_5p_end', 'peak2_score', 'peak2_signal_value', 'peak2_pval', 'peak2_qval', 'peak2_point_source', 'peak2_start', 'peak2_end', 'distance', 'distance_signed']
    df[numeric_cols] = df[numeric_cols].fillna(value = 0)
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)
    
    return df


def filter_bipeak_dist(row, row_ind, df, dist):
    
    peak_dist = row.distance
    is_bipeak = str(peak_dist) < str(dist)
    
    return is_bipeak


def filter_df(df, dist):
    
    bipeaks = []
    unipeaks = []
    
    for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Filter df"):
        row = row_vals[1]
        if pd.isnull(row.peak1_id) or pd.isnull(row.peak2_id):
            row.distance = row.bipeak_end - row.bipeak_start
        if filter_bipeak_dist(row, row_ind, df, dist):
            bipeaks.append(row)
        else:
            unipeaks.append(row)
            
    return pd.DataFrame(bipeaks), pd.DataFrame(unipeaks)


def unique_bipeaks(df, dedup_cols):
    
    dedup = df[0].drop_duplicates(subset=dedup_cols)
    dup1 = dedup.duplicated(subset=['peak1_id'], keep=False)
    peak1 = dedup[dup1]
    dup2 = dedup.duplicated(subset=['peak2_id'], keep=False)
    peak2 = dedup[dup2]
    bipeaks = dedup[~(dup1 | dup2)]
    
    prev_row = None
    for row_ind, row_vals in tqdm(enumerate(peak1.iterrows()), total=len(peak1), desc="Unique bipeaks 1"):
        row = row_vals[1]
        if prev_row is not None and row.peak1_id == prev_row.peak1_id:
            if row.distance < prev_row.distance:
                bipeaks.append(row)
                # bipeaks = pd.concat([bipeaks, row]).reset_index(drop=True)
            else:
                bipeaks.append(prev_row)
                # bipeaks = pd.concat([bipeaks, prev_row]).reset_index(drop=True)
        prev_row = row
    
    prev_row = None
    for row_ind, row_vals in tqdm(enumerate(peak2.iterrows()), total=len(peak2), desc="Unique bipeaks 2"):
        row = row_vals[1]
        if prev_row is not None and row.peak2_id == prev_row.peak2_id:
            if row.distance < prev_row.distance:
                bipeaks.append(row)
                # bipeaks = pd.concat([bipeaks, row]).reset_index(drop=True)
            else:
                bipeaks.append(prev_row)
                # bipeaks = pd.concat([bipeaks, prev_row]).reset_index(drop=True)
        prev_row = row
        
    print("\nNo peak 1 duplicates remaining:", not any(bipeaks.duplicated(subset = 'peak1_id', keep = False)))
    print("No peak 2 duplicates remaining:", not any(bipeaks.duplicated(subset = 'peak2_id', keep = False)))
    
    return pd.DataFrame(bipeaks)


def extract_unique_peaks(df, peak_list, peak_id):
    
    peak_ind = 0
    peak_df = []
    
    for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df), desc="Extract unique peaks"):
        row = row_vals[1]
        if peak_ind >= len(peak_list):
            break
        if row[peak_id] == peak_list[peak_ind]:
            peak_df.append(row)
            peak_ind += 1
            
    return pd.DataFrame(peak_df)


def generate_unipeak_ids(unipeaks):
    
    cell_line = unipeaks.peak_id[0].split(".", 1)[0]
    unipeak_ids = []
    
    for i in tqdm(range(1, len(unipeaks)+1), total=len(unipeaks), desc="Generate ids"):
        unipeak_id = cell_line + "." + "unipeak_" + str(i)
        unipeak_ids.append(unipeak_id)
        
    return unipeak_ids


def unique_unipeaks(df, bipeaks, dedup_cols):
    
    dedup = df[1].drop_duplicates(subset=dedup_cols)
    
    #create sets of bipeak and unipeak peak1 and peak2 values and get intersection with unipeak - bipeak
    bipeaks1 = set(bipeaks.peak1_id)
    bipeaks2 = set(bipeaks.peak2_id)
    unipeaks1 = {x for x in set(dedup.peak1_id) if pd.notna(x)}
    unipeaks2 = {x for x in set(dedup.peak2_id) if pd.notna(x)}
    peak1 = unipeaks1 - bipeaks1
    peak2 = unipeaks2 - bipeaks2
    
    #extract one row with each unique peak_id from original dataframe
    peak1_list = sorted(peak1)
    peak2_list = sorted(peak2)
    unique1 = extract_unique_peaks(dedup.sort_values(by = 'peak1_id'), peak1_list, 'peak1_id')
    unique2 = extract_unique_peaks(dedup.sort_values(by = 'peak2_id'), peak2_list, 'peak2_id')
    
    #reformat dataframe
    unipeaks1 = unique1[['peak1_chr', 'peak1_start', 'peak1_end', 'peak1_score', 'peak1_strand', 'peak1_id', 'peak1_5p_start', 'peak1_5p_end', 'peak1_signal_value', 'peak1_pval', 'peak1_qval', 'peak1_point_source']]
    unipeaks2 = unique2[['peak2_chr', 'peak2_start', 'peak2_end', 'peak2_score', 'peak2_strand', 'peak2_id', 'peak2_5p_start', 'peak2_5p_end', 'peak2_signal_value', 'peak2_pval', 'peak2_qval', 'peak2_point_source']]
    unipeaks1.columns = ['unipeak_chr', 'unipeak_start', 'unipeak_end', 'unipeak_score', 'unipeak_strand', 'peak_id', '5p_start', '5p_end', 'signal_value', 'pval', 'qval', 'point_source']
    unipeaks2.columns = ['unipeak_chr', 'unipeak_start', 'unipeak_end', 'unipeak_score', 'unipeak_strand', 'peak_id', '5p_start', '5p_end', 'signal_value', 'pval', 'qval', 'point_source']
    unipeaks = pd.concat([unipeaks1, unipeaks2], ignore_index=True)
    unipeaks = unipeaks.sort_values(by = ['unipeak_chr', 'unipeak_start'])
    unipeaks.reset_index(drop=True, inplace=True)
    unipeak_ids = generate_unipeak_ids(unipeaks)
    class_name = ['single'] * (len(unipeaks))
    unipeaks.insert(3, "unipeak_id", pd.Series(unipeak_ids))
    unipeaks['unipeak_class'] = class_name
    
    #assess peak_id uniqueness
    print("\nNo unipeak duplicates remaining:", not any(unipeaks.duplicated(subset = 'peak_id', keep = False)))
    
    return unipeaks.rename(columns={'unipeak_chr': '#unipeak_chr'})


def create_final_df(bipeaks, unipeaks):
    
    #subset and reformat bipeak file
    bp = bipeaks["peak_id"] = bipeaks["peak1_id"] + "," + bipeaks["peak2_id"]
    bp = bipeaks[['bipeak_chr', 'bipeak_start', 'bipeak_end', 'bipeak_id', 'bipeak_score', 'bipeak_strand', 'peak_id', 'bipeak_midpoint', 'distance', 'bipeak_class']]
    bp.columns = ['chr', 'start', 'end', 'class_id', 'score', 'strand', 'peak_id', 'bipeak_midpoint', 'distance', 'class']
    
    #subset and reformat unipeak file
    up = unipeaks[['unipeak_chr', 'unipeak_start', 'unipeak_end', 'unipeak_id', 'unipeak_score', 'unipeak_strand', 'peak_id', 'unipeak_class']]
    up_empty = ['.'] * (len(up))
    up.insert(7, "bipeak_midpoint", pd.Series(up_empty))
    up.insert(8, "distance", pd.Series(up_empty))
    up.columns = ['chr', 'start', 'end', 'class_id', 'score', 'strand', 'peak_id', 'bipeak_midpoint', 'distance', 'class']
    
    #concatenate dataframes and sort
    allpeaks0 = pd.concat([bp, up])
    allpeaks = allpeaks0.sort_values(by = ['chr', 'start']).rename(columns={'chr': '#chr'})
    
    return allpeaks


def main(bed_filepath, dist):
    
    header = ['bipeak_chr', 'bipeak_start', 'bipeak_end', 'bipeak_id', 'bipeak_score', 'bipeak_strand', 'bipeak_midpoint', 'peak1_chr', 'peak1_5p_start', 'peak1_5p_end', 'peak1_id', 'peak1_score', 'peak1_strand', 'peak1_signal_value', 'peak1_pval', 'peak1_qval', 'peak1_point_source', 'peak1_start', 'peak1_end', 'peak2_chr', 'peak2_5p_start', 'peak2_5p_end', 'peak2_id', 'peak2_score', 'peak2_strand', 'peak2_signal_value', 'peak2_pval', 'peak2_qval', 'peak2_point_source', 'peak2_start', 'peak2_end', 'distance', 'distance_signed', 'bipeak_class']
    dedup_cols = header[:3] + header[4:]
    
    print("\nDownloading dataframe")
    df = load_from_bed(bed_filepath, header)
    
    print("\nFiltering dataframe")
    filtered_df = filter_df(df, dist)
    
    print("\nGenerating bipeak file")
    bipeaks0 = unique_bipeaks(filtered_df, dedup_cols)
    bipeaks = bipeaks0.rename(columns={'bipeak_chr': '#bipeak_chr'})
    
    print("\nGenerating unipeak file")
    unipeaks0 = unique_unipeaks(filtered_df, bipeaks, dedup_cols) #generate_ids not working? No seems to be... must be merge that is wrong?
    unipeaks = unipeaks0.rename(columns={'unipeak_chr': '#unipeak_chr'})
    
    print("\nGenerating merged all_peak file")
    all_peaks = create_final_df(bipeaks, unipeaks)

    #Output files
    cell_name = ".".join(bed_filepath.split(".")[:-2])
    bipeak_output = cell_name+".unique_bipeaks_"+str(dist)+".bed"
    unipeak_output = cell_name+".unique_unipeaks_"+str(dist)+".bed"
    allpeak_output = cell_name+".unique_peaks_"+str(dist)+".bed"
    
    delete_if_exists(bipeak_output)
    delete_if_exists(unipeak_output)
    delete_if_exists(allpeak_output)
    
    bipeaks.to_csv(bipeak_output, sep='\t', header=True, index=False)
    unipeaks.to_csv(unipeak_output, sep='\t', header=True, index=False)
    all_peaks.to_csv(allpeak_output, sep='\t', header=True, index=False)


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('file',
                        help='Path to file.')
    parser.add_argument('--dist', default=1000,
                        help='Distance threshold between two paired peaks. If not specified, defaults to 1kb.')
    
    args = parser.parse_args()
   
    main(args.file, args.dist)
	

