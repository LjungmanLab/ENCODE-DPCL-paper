# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 16:07:48 2022

@author: abmcs
"""

import argparse
import uuid
import os
import pandas as pd
import re
from collections import defaultdict


def filter_list(datalist):
    
    return [val for val in datalist if re.search(r'-all.peak_overlaps.bed', val)]


def delete_if_exists(filepath):
    
    if os.path.exists(filepath):
        os.remove(filepath)
        return True
    
    return False


def load_files_to_dict(bed_filepath_list):
    
    header = ['chr1', 'start1', 'end1', 'class_id1', 'score1', 'strand1', 'peak_id1', 'bipeak_midpoint1', 'distance1', 'class1', 'num_b_ovlps', 'chr2', 'start2', 'end2', 'class_id2', 'score2', 'strand2', 'peak_id2', 'bipeak_midpoint2', 'distance2', 'class2', 'bp_ovlp']
    
    dfs = {}
    for filepath in bed_filepath_list:
        df = pd.read_csv(filepath, sep='\t', names=header, engine='python')
        cell_line = df.iloc[0,3].split(".")[0]
        dfs[cell_line] = df
    
    return dfs


def get_peak_id_batches(df):
    """
    Parameters
    ----------
    df : pandas dataframe containing results from bedtools intersection of cell_line peaks
         vs all other peaks (i.e. 'A673-all.peak_overlaps.bed')

    Returns
    -------
    A list of lists of peak ids (just strings) associated with the same group
    
    Ex: 
        [["GM12878.bipeak_129", "HCT116.unipeak_106", ...], ["Caco2.unipeak_70", ...], ...]
    """

    unique_peak_ids = set(df.class_id1)
    
    id_batches = []
    for peak_id in unique_peak_ids:
        sub_df = df[df['class_id1'] == peak_id]
        if len(sub_df) == 1 and (sub_df['class_id2'] == '.').bool():
            peak_id_list = [peak_id]
        else:
            peak_id_list = sub_df.class_id2.values.tolist()
            peak_id_list.append(peak_id)
        
        id_batches.append(peak_id_list)
    
    return id_batches


def merge_groups(peakid_to_group, group_to_peakid, peak_group_ids):
    '''
    Parameters
    ----------
    peakid_to_group : dictionary where keys = peak_ids and values = group uuids
    group_to_peakid : dictionary where values = group uuids and keys = list of overlapping peak_ids
    peak_group_ids : set containing all group uuids associated with peak ids being looped through

    Returns
    -------
    Both dictionaries detailed above (with old keys/values removed) and a list of peak_ids in the 
    group that will eventually be added as a value to group uuid key in group_to_peakid
    '''
    
    # Get all peak ids from each group in peak_group_ids
    # Add all peaks ids to a new set() (peak_ids)
    # For each group id, remove it from group_to_peakid
    peak_ids = set()    
    for group_id in peak_group_ids:
        id_list = group_to_peakid.get(group_id)
        peak_ids.update(id_list)
        group_to_peakid.pop(group_id)
    
    # For each peak id, remove it from peakid_to_group
    for peak_id in list(peak_ids):
        peakid_to_group.pop(peak_id)
    
    return peakid_to_group, group_to_peakid, list(peak_ids)


def compile_peak_id_groups(df_dict):
    '''
    Parameters
    ----------
    df_dict : dictionary where keys = cell line and values = dataframes containing results from 
              bedtools intersection of cell_line peaks vs all other peaks (i.e. 'A673-all.peak_overlaps.bed') 

    Returns
    -------
    peakid_to_group, a dictionary where keys = peak_ids and values = group uuids and group_to_peakid, a 
    dictionary where values = group uuids and keys = list of overlapping peak_ids
    '''
    
    peakid_to_group = {}
    group_to_peakid = {}
    for df in df_dict.values():
        # get rows with same peak_ids in one file and get list of all overlapping peak_ids
        for peak_ids in get_peak_id_batches(df):
            # check in peak_ids already exist in peakid_to_group
            peak_group_ids = set()
            for peak_id in peak_ids:
                if peak_id in peakid_to_group:
                    group_id = peakid_to_group[peak_id]
                    peak_group_ids.add(group_id)
            # if yes, get group_ids and either remove group_ids from group_to_peakid and reassign (if mulitple)
            if len(peak_group_ids) > 1:
                peakid_to_group, group_to_peakid, peak_ids = merge_groups(peakid_to_group, group_to_peakid, peak_group_ids)
                group_id = str(uuid.uuid4())
                peak_group_ids = set()
            # or, get group_id and use it for new peak_id
            if len(peak_group_ids) == 1:
                group_id = peak_group_ids.pop()
            # if no, generate new group_id
            else:
                group_id = str(uuid.uuid4())
            
            # then, assign all peak_ids to group_id and add to group_to_peakid and peakid_to_group
            existing_peak_ids = group_to_peakid.get(group_id, [])
            group_to_peakid[group_id] = list(set(existing_peak_ids+peak_ids))
            for peak_id in peak_ids:
                peakid_to_group[peak_id] = group_id
        
    return peakid_to_group, group_to_peakid



######### end df parsing code #########
########### start eval code ###########



def find_cell_matches(peak_id_list):
    '''
    Parameters
    ----------
    peak_id_list : list of peak ids

    Returns
    -------
    A dictionary where keys = cell lines and values = peak_id lists (used in get_groups_to_eval())
    '''
    
    cell_dict = {}
    for peak_id in peak_id_list:
        cell_line = peak_id.split(".")[0]
        if cell_line in cell_dict:
            cell_dict[cell_line]+=[peak_id]
        else:
            cell_dict[cell_line] = [peak_id]
    
    return cell_dict


def get_groups_to_eval(eval_groups, cell_dict, group_id):
    '''
    Parameters
    ----------
    eval_groups : empty list
    cell_dict : dictionary output from find_cell_matches()
    group_id : group_id of peak group currently being evaluated

    Returns
    -------
    A list containing group_ids that need to be parsed
    '''
    
    # if more than 1 peak id in cell_dict group requires eval and group_id is added to eval groups
    for cell in cell_dict:
        if len(cell_dict[cell]) > 1:
            eval_groups.append(group_id)            
            break
    
    return eval_groups


def get_df_lists_to_eval(fixed_group_df_list):
    '''
    Parameters
    ----------
    fixed_group_df_list : list of dataframes that underwent a primary evaluation (primary_eval_function)

    Returns
    -------
    duplicate_dfs, a list of dataframes that still contain multiple peak_ids from the same cell line and
    finalized_dfs, a list of dataframes that no longer contain multiple peak_ids from the same cell line
    '''
    
    duplicate_dfs = []
    finalized_dfs = []
    for fixed_group_df in fixed_group_df_list:
        # verify that there are no more duplicates in these dfs
        fixed_pkid_list = list(fixed_group_df.class_id1)
        fixed_cell_dict = find_cell_matches(fixed_pkid_list)
    
        # if more than 1 peak id in cell_dict group requires eval and group_id is added to eval groups
        is_duplicate = False    
        for cell in fixed_cell_dict:
            if len(fixed_cell_dict[cell]) > 1:
                is_duplicate = True
        if is_duplicate:
            duplicate_dfs.append(fixed_group_df)
        else:
            finalized_dfs.append(fixed_group_df)
    
    return duplicate_dfs, finalized_dfs


def compile_group_dfs(group_dict, df_dict):
    '''
    Parameters
    ----------
    group_dict : dictionary where keys = group_ids and values = peak_ids
    df_dict : dictionary where keys = cell line and values = dataframes containing results from 
              bedtools intersection of cell_line peaks vs all other peaks (i.e. 'A673-all.peak_overlaps.bed') 

    Returns
    -------
    A dictionary where keys are group_ids and values are dataframes containing rows that coorespond to the 
    peak_ids in the group and columns 1-10 of the original df (from df_dict)
    '''
    group_df_dict = {}
    for group in group_dict:
        
        group_df = []
        for peak in group_dict[group]:
            cell_line = peak.split(".")[0]
            df = df_dict[cell_line]
            row = df[df.class_id1 == peak].iloc[0,0:10]
            group_df.append(row)
        
        pd_group_df = pd.DataFrame(group_df)
        pd_group_df = pd_group_df.reset_index(drop = True)
        
        group_df_dict[group] = pd_group_df

    return group_df_dict


def primary_eval_function(df):
    '''
    *** BEWARE: Recursive Function. Loop carefully! ***

    Parameters
    ----------
    df : pandas dataframe

    Returns
    -------
    A list of dataframes that have been evaluated and subset to contain peaks that are no more than 600bp
    from the midpoint of all peaks in the group (some exceptions exist).

    '''
    # end recursive function if input df is empty
    if len(df) == 0:
        return []
    
    # create new column in df with coordinates to be aligned (eval_col)
    eval_col = []
    for row_ind, row_vals in df.iterrows():
        row = row_vals
        if row.bipeak_midpoint1 != '.':
            eval_col.append(int(row.bipeak_midpoint1))
        if row.strand1 == '+':
            eval_col.append(int(row.start1))
        if row.strand1 == '-':
            eval_col.append(int(row.end1))
    
    # add eval_col 
    df.loc[:,"eval_coordinate"] = eval_col
    
    # find the distance of that coordinate from the median of the column
    df.loc[:,"dist_from_median"] = df.loc[:,"eval_coordinate"].median() - df.loc[:,"eval_coordinate"]
    
    # if coordinate is more than 600bp away from median remove from group
    fixed_group_df = df.drop(df[abs(df.dist_from_median) > 600].index)
    
    # get the other peak_ids that need to be removed from group
    new_eval_df = df.loc[abs(df.dist_from_median) > 600]


    if len(new_eval_df) == 2 and "GM12878.unipeak_42655" in new_eval_df.class_id1 and "A673.unipeak_25284" in new_eval_df.class_id1:
        import pdb; pdb.set_erace()
    ### check for exceptions:
    # if two peaks left in df that are not within 600bp of each other, return each as individual df
    if len(new_eval_df) == len(df) == 2:
        return [df.iloc[[0]], df.iloc[[1]]]
    # if more than two peaks left that are not within 600bp of median, check strand
    if len(new_eval_df) == len(df) > 2:
            #if all unipeaks and both strands exist, split on strand
        try:
            if '.' not in list(df['strand1']) and not (df['strand1'] == df['strand1'].iloc[0]).all():
                plus_df = df.loc[df.strand1 == '+'].copy()
                minus_df = df.loc[df.strand1 == '-'].copy()
                recursive_fixed_group_dfs_1 = primary_eval_function(plus_df)
                recursive_fixed_group_dfs_2 = primary_eval_function(minus_df)
                recursive_fixed_group_dfs = recursive_fixed_group_dfs_1 + recursive_fixed_group_dfs_2
                return recursive_fixed_group_dfs
            #if bipeak or all strands are same, return df
            else:
                return [new_eval_df]
        except KeyError:
            import pdb; pdb.set_trace()
    
    return [fixed_group_df] + primary_eval_function(new_eval_df)


def secondary_eval_function(duplicate_dfs, finalized_dfs):
    # Rules for secondary eval:
            # If no bipeak, sep +/- peaks
            # If 1 bi + 1 uni for cell line, keep bi and put uni in new group
            # Else leave as group
    for dup_df in duplicate_dfs:
                # split on "." and "_" and add cell line and peak type info to dictionary
                peakid_breakdown = {}
                cell_line_type_matches = defaultdict(list)
                all_unipeaks = True
                for row_ind, row_vals in dup_df.iterrows():
                    row = row_vals
                    cell_line = re.split('\.|_', row.class_id1)[0]
                    peak_type = re.split('\.|_', row.class_id1)[1]
                    peakid_breakdown[str(row.class_id1)] = (cell_line, peak_type)
                    # set to false if any bipeaks in the df
                    if peak_type == 'bipeak':
                        all_unipeaks = False
                    cell_line_type_matches[cell_line].append((str(row.class_id1), peak_type))
            
                # check if any bipeaks and if no get df and split on strand col
                if all_unipeaks:
                    plus_df = dup_df[dup_df.strand1 == '+']
                    minus_df = dup_df[dup_df.strand1 == '-']
                    finalized_dfs.append(plus_df)
                    finalized_dfs.append(minus_df)
                # check if there is a bipeak and unipeak for one cell line
                else:
                    matched_unipeak_ids = []
                    for cell_line, peak_tuples in cell_line_type_matches.items():
                        types = [p[1] for p in peak_tuples]
                        if "unipeak" in types and "bipeak" in types:
                            for p in peak_tuples:
                                peak_id, peak_type = p
                                if peak_type == "unipeak":
                                    matched_unipeak_ids.append(peak_id)
                    # if no, keep original df and add to finalized_dfs
                    if not matched_unipeak_ids:
                        finalized_dfs.append(dup_df)
                    # if yes, pull out unipeak and anything with midpoint withing 600bp of it
                    # then reevaluate all resulting dfs
                    else:
                        for unipeak in matched_unipeak_ids:
                            eval_coord = dup_df.query('class_id1 == @unipeak')['eval_coordinate']
                            unipeak_df = dup_df.loc[abs(dup_df.loc[:,"eval_coordinate"] - int(eval_coord)) < 600]
                            bipeak_df = dup_df.loc[abs(dup_df.loc[:,"eval_coordinate"] - int(eval_coord)) > 600]
                            if len(unipeak_df) > 0:
                                finalized_dfs.append(unipeak_df)
                            if len(bipeak_df) > 0:
                                finalized_dfs.append(bipeak_df)
                        
    return finalized_dfs


def eval_testing_function(evaluated_dfs):
    '''
    Parameters
    ----------
    evaluated_dfs : TYPE
        DESCRIPTION.

    Returns
    -------
    output_dfs : TYPE
        DESCRIPTION.
    '''
    
    output_dfs = []
    for eval_df in evaluated_dfs:
        ### TEST again with primary and secondary
        actual_final_df_list = primary_eval_function(eval_df)
        duplicate_dfs, finalized_dfs = get_df_lists_to_eval(actual_final_df_list)
        if duplicate_dfs:
            finalized_dfs = secondary_eval_function(duplicate_dfs, finalized_dfs)
            
        output_dfs += finalized_dfs
        
    return output_dfs

def are_equal(df_list1, df_list2):
    if len(df_list1) != len(df_list2):
        return False
    
    for df1, df2 in zip(df_list1, df_list2):
        if df1.shape != df2.shape or not all(df1 == df2):
            return False
        
    return True


def evaluate_peaks(eval_dfs, peakid_to_group, group_to_peakid):
    '''
    Parameters
    ----------
    eval_dfs : dictionary where keys = group ids and values = dataframes to be evaluated
    peakid_to_group : dictionary where keys = peak_ids and values = group uuids
    group_to_peakid : dictionary where values = group uuids and keys = list of overlapping peak_ids

    Returns
    -------
    The updated peakid_to_group and group_to_peakid dictionaries
    '''
    
    for df_id in eval_dfs:
        df = eval_dfs[df_id]
        # primary evaluation of groups based on proximity to median coordinate
        fixed_group_df_list = primary_eval_function(df)
        # check groups after first eval and get a list of fixed dfs and dfs that need to be reevaluated
        duplicate_dfs, evaluated_dfs = get_df_lists_to_eval(fixed_group_df_list)
        
        # if df that need further eval in duplicate_dfs, perform secondary eval on these
        if duplicate_dfs:
            evaluated_dfs = secondary_eval_function(duplicate_dfs, evaluated_dfs)
        
        # test evaluated dfs to ensure that all exceptions are taken care of
        finalized_dfs = eval_testing_function(evaluated_dfs)
        while not are_equal(finalized_dfs, evaluated_dfs):
            evaluated_dfs = finalized_dfs
            finalized_dfs = eval_testing_function(evaluated_dfs)
        
        # update group_to_peakid and peakid_to_group
        for finalized_df in finalized_dfs:
            if len(finalized_df) != 0:
                peak_id_list = list(finalized_df.class_id1)
                final_group_id = str(uuid.uuid4())
                group_to_peakid[final_group_id] = peak_id_list
                for peak_id in peak_id_list:
                    peakid_to_group[peak_id] = final_group_id

    return peakid_to_group, group_to_peakid


def fix_peak_groups(peakid_to_group, group_to_peakid, df_dict):
    '''
    Parameters
    ----------
    peakid_to_group : dictionary where keys = peak_ids and values = group uuids
    group_to_peakid : dictionary where values = group uuids and keys = list of overlapping peak_ids
    df_dict : dictionary where keys = cell line and values = dataframes containing results from 
              bedtools intersection of cell_line peaks vs all other peaks (i.e. 'A673-all.peak_overlaps.bed') 

    Returns
    -------
    The updated peakid_to_group (new_peakid_to_group) and group_to_peakid dictionaries (new_group_to_peakid)
    '''
    
    # determine if group in group_to_peakid has multiple peakids from a single cell line and collect those group_ids in a list
    eval_groups = []
    for group_id in group_to_peakid:
        # if only one peak_id in group, do not check
        if len(group_to_peakid[group_id]) == 1:
            continue
        # get peak_ids in a list
        peak_id_list = group_to_peakid[str(group_id)]
        # create dictionary to find multiple peaks from one cell line in group
        cell_dict = find_cell_matches(peak_id_list)
        # add group to list if mult peaks from cell line and needs to be evaluated
        eval_groups = get_groups_to_eval(eval_groups, cell_dict, group_id)
    
    # create dictionary of groups to peakids for eval_groups
    eval_dict = {key: group_to_peakid[key] for key in eval_groups}
    
    # remove group_id from group_to_peakid and remove peak_ids from peakid_to_group
    for eval_group_id, peak_list in eval_dict.items():
        group_to_peakid.pop(eval_group_id)
        [peakid_to_group.pop(peak) for peak in peak_list]
    
    # create dictionary of dfs where each row has the information for each peak_id in the eval_group
    eval_dfs = compile_group_dfs(eval_dict, df_dict)
    
    # remove duplicate peaks in each group and reassign groups to
    new_peakid_to_group, new_group_to_peakid = evaluate_peaks(eval_dfs, peakid_to_group, group_to_peakid)
    
    return new_peakid_to_group, new_group_to_peakid


def list2Str(lst):
    if type(lst) is list: # apply conversion to list columns
        return",".join(lst)
    else:
        return lst
    

def assemble_data_matrix(group_to_peakid, df_dict):
    '''
    Parameters
    ----------
    group_to_peakid : dictionary where values = group uuids and keys = list of overlapping peak_ids
    df_dict : dictionary where keys = cell line and values = dataframes containing results from 
              bedtools intersection of cell_line peaks vs all other peaks (i.e. 'A673-all.peak_overlaps.bed') 

    Returns
    -------
    A dataframe where rows are newly defined peak regions and columns are cell lines, and cell values
    are the peak ids for the cell line that fall in the new peak region
    '''
    
    # make dictionary of dfs for all groups
    group_df_dict = compile_group_dfs(group_to_peakid, df_dict)
    
    # generate new group coordinates and create dictionary linking group_id and coordinate
    group_to_coord = {}
    for group_id, df in group_df_dict.items():
        chrom = df.iloc[0,0]
        start = df['start1'].min()
        end = df['end1'].max()
        new_id = chrom+":"+str(start)+"-"+str(end)
        group_to_coord[group_id] = new_id
    
    # create dictionary of coordinates linked to peak id lists
    coord_to_peakid = dict((group_to_coord[key], value) for (key, value) in group_to_peakid.items())
    
    # create nested dictionary as framework for dataframe
    df_framework_dict = {}
    for coord, peak_id_list in coord_to_peakid.items():
        peak_id_dict = {}
        for peak_id in peak_id_list:
            cell_line = peak_id.split(".")[0]
            if cell_line in peak_id_dict:
                peak_id_dict[cell_line] += [peak_id]
            else:
                peak_id_dict[cell_line] = [peak_id]
        df_framework_dict[coord] = peak_id_dict
    
    # create dataframe from dictionary
    pd_df = pd.DataFrame.from_dict(df_framework_dict, orient='index')
    
    # convert lists to strings
    final_df = pd_df.apply(lambda x: [list2Str(i) for i in x])
    
    # move index to column
    final_df.reset_index(inplace = True)
    final_df = final_df.rename(columns = {'index':'coordinate', 'A673UV':'A673', 'Caco2UV':'Caco2', 
                                          'Calu3UV':'Calu3', 'GM12878UV':'GM12878', 'HCT116UV':'HCT116', 
                                          'HepG2UV':'HepG2', 'IMR90UV':'IMR90', 'K562UV':'K562',
                                          'MCF10AUV':'MCF10A', 'MCF7UV':'MCF7', 'OCILY7UV':'OCILY7',
                                          'panc1UV':'panc1','PC3UV':'PC3', 'PC9UV':'PC9'})
    final_df = final_df.reindex(['coordinate','A673','Caco2','Calu3','GM12878','HCT116','HepG2','IMR90','K562','MCF10A','MCF7','OCILY7','panc1','PC3','PC9'], axis = 1)
    
    return final_df


def create_final_bed(df):
    
    df.insert(1, 'strand', '.')
    df.insert(1, 'score', '.')
    
    chrom = df['coordinate'].str.split(':').str[0]
    coord = df['coordinate'].str.split(':').str[1]
    start = coord.str.split('-').str[0]
    end = coord.str.split('-').str[1]
    
    df.insert(0, 'end', pd.to_numeric(end))
    df.insert(0, 'start', pd.to_numeric(start))
    df.insert(0, 'chr', chrom)
    
    final_bed = df.sort_values(["chr", "start"])
    final_bed = final_bed.reset_index(drop=True)
    
    return final_bed


def main(bed_filepath_list):
    
    df_dict = load_files_to_dict(bed_filepath_list)
    peak_ids, group_ids = compile_peak_id_groups(df_dict)
    
    #evaluate peak id lists in each group and resolve instances of multiple peaks from same file
    final_peak_ids, final_group_ids = fix_peak_groups(peak_ids.copy(), group_ids.copy(), df_dict)
    
    #then make dataframe
    final_df = assemble_data_matrix(final_group_ids, df_dict)
    
    #make bed file format
    final_bed = create_final_bed(final_df)
    
    #write files
    tsv_file_name = "BruUV_eRNA_16_cell_lines.tsv"
    delete_if_exists(tsv_file_name)
    final_bed.to_csv(tsv_file_name, sep='\t', na_rep="NA", header=True, index=False)

    bed_file_name = "BruUV_eRNA_16_cell_lines.bed"
    delete_if_exists(bed_file_name)
    final_bed.to_csv(bed_file_name, sep='\t', na_rep="NA", header=True, index=False)


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('directory',
                        help='Path to file directory.')
    
    args = parser.parse_args()
    
    filepath_list = [os.path.abspath(os.path.join(args.directory, p)) for p in os.listdir(args.directory)]
    bed_filepath_list = filter_list(filepath_list)
    
    main(bed_filepath_list)
    