# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:28:01 2022

Description
-----------
This script defines the following functions which are used to... :
    
    load_from_bed: loads in bed file as pandas compatible dataframe and adds header
    delete_if_exists: deletes file if already exists in working directory
    subset_rpkm: gets rpkm column from basic rpkm file to merge with RT_segment df
    calc_gene_bp: calculates gene lengths to add to RT_segment df
    redefine_rt_seg: redefines segment end coordinate as gene tes
    calc_tes_seg_diff: recalculates RT_length (distance from tes to segment end)
    reformat_df: adds new columns to and reformats RT_segment df
    filter_genes: filters RT_segments according to segment score, gene length, and gene rpkm
    assign_rt_segments: reassigns duplicate RT_segments (overlap multiple genes) to the most upstream gene TES
    
@author: abmcs
"""

# load data as pandas dataframe

import argparse
import os
from tqdm import tqdm
import pandas as pd


def load_from_bed(filepath):
    
    '''
    loads in fused_dec file (sample.segments.fused_dec_simple.tes.intersect_test.bed) output
    from intersect with tes_annot (see bash rt_workingscript) and adds column headers
    '''
    
    header = ['gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'min_segment_score', 'num_fused_seg', 'intersect_bp']
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


def subset_rpkm(counts_df, rpkm_df, sample):
    
    '''
    splits gene column into ID and name and subsets rpkm dataframe to contain gene_ID,
    gene_name, and sample (gene_rpkm) column for merging with read-through df
    '''
    
    counts_df = counts_df.rename(columns={'geneName/geneID': 'gene', 'featureLength': 'length'})
    rpkm_df = rpkm_df.rename(columns={'geneName/geneID': 'gene', 'featureLength': 'length'})
    
    count_col = 'gene_counts_' + sample
    
    counts_df[['gene_name', 'gene_ID']] = counts_df['gene'].str.rsplit('/', n=1, expand=True)
    counts_df['density'] = counts_df[count_col] / counts_df['length']
    counts_sub = counts_df[['gene_ID', 'gene_name', count_col, 'density']]
    
    rpkm_col = 'gene_RPKM_' + sample
    
    rpkm_df[['gene_name', 'gene_ID']] = rpkm_df['gene'].str.rsplit('/', n=1, expand=True)
    rpkm_sub = rpkm_df[['gene_ID', 'gene_name', rpkm_col]]
    
    merge_df = counts_sub.merge(rpkm_sub, how='outer', on=['gene_ID', 'gene_name'])
    
    out_df = merge_df.rename(columns={rpkm_col:'gene_rpkm', count_col:'gene_counts', 'density':'gene_density'})
    
    return out_df


def calc_gene_bp(row):
    
    '''
    calculates the gene length and reports an integer that will be added as a
    column in reformatted df
    '''
    
    gene_length = row.gene_end - row.gene_start
    
    return gene_length


def redefine_rt_seg(row, new_seg_start, new_seg_end):
    
    '''
    based on strand, redefines the RT_segment start coordinate as the gene tes
    and the RT_segment end coordinate as the end of the segment
    '''
    
    if row.segment_strand == '-':
        row[new_seg_start] = row.segment_start
        row[new_seg_end] = row.tes_start
    else:
        row[new_seg_start] = row.tes_end
        row[new_seg_end] = row.segment_end
    
    return


def calc_tes_seg_diff(row):
    
    '''
    calculates RT_length (distance from tes to segment end) and reports an
    integer that will be added as a column in reformatted df
    '''
    
    if row.segment_strand == '-':
        ts_diff = row.tes_start - row.segment_start
    else:
        ts_diff = row.segment_end - row.tes_end
    
    return int(ts_diff)


def reformat_df(df, rpkm_df):
    
    '''
    Parameters
    ----------
    df: takes fused_dec output that has been intersected with your tes annotation
        (suffix = .segments.fused_dec_simple.tes.intersect_test.bed)
    rpkm_df: takes gene rpkm output from make_counts_matrix_PEdata.sh script
        containing standard bed columns plus sample RPKM columns (expects column header
        to == sample_name, and expects gene column to contain ensembl_id/gene_name)
    
    Description
    -----------
    This function takes read-through segments intersecting TESs and a gene rpkm matrix
    and reformats the read-through matrix such that:
        - genes/tes that do not intersect a segment are removed
        - RT_segment coordinates are redefined as the tes -> start of segment
        - new columns are generated to describe RT_segments
        - gene rpkm and gene_length columns are added
    
    Variables
    ---------
    No variables require changing but be aware of df formatting/column expectations/requirements
    '''
    
    new_df = []
    tes_seg_diff = 'RT_length'
    new_seg_start = 'RT_start'
    new_seg_end = 'RT_end'
    new_seg_id = 'RT_ID'
    new_seg_feat = 'RT_feat_type'
    new_seg_score = 'RT_score'
    new_seg_strand = 'RT_strand'
    new_seg_chr = 'RT_chr'
    gene_len = 'gene_length'
    
    for row_ind, row_vals in tqdm(enumerate(df.iterrows()), total=len(df)):
        row = row_vals[1]
        seg_chr = row.segment_chr
        if seg_chr != '.':
            row[tes_seg_diff] = calc_tes_seg_diff(row)
            row[gene_len] = calc_gene_bp(row)
            redefine_rt_seg(row, new_seg_start, new_seg_end)
            new_df.append(row)
            
    pd_df = pd.DataFrame(new_df)
    pd_df[new_seg_chr] = pd_df['segment_chr']
    pd_df[new_seg_id] = pd_df['RT_chr']+":"+pd_df['RT_start'].astype(str)+"-"+pd_df['RT_end'].astype(str)
    pd_df[new_seg_score] = pd_df['segment_score']
    pd_df[new_seg_strand] = pd_df['segment_strand']
    pd_df[new_seg_feat] = 'RT_segment'
    
    rpkm_df = pd_df.merge(rpkm_df, how='left', on=['gene_ID', 'gene_name'])
    
    reform_df = rpkm_df[['RT_chr', 'RT_start', 'RT_end', 'RT_ID', 'RT_score', 'RT_strand', 'RT_feat_type', 'RT_length', 'gene_chr', 'tes_start', 'tes_end', 'gene_ID', 'gene_score', 'gene_strand', 'gene_name', 'gene_feat_type', 'gene_start', 'gene_end', 'gene_length', 'gene_counts', 'gene_density', 'gene_rpkm', 'segment_chr', 'segment_start', 'segment_end', 'segment_ID', 'segment_score', 'segment_strand', 'segment_feat_type', 'segment_length', 'fused_segment_ids', 'min_segment_score', 'num_fused_seg', 'intersect_bp']]
    
    numeric_cols = ['RT_start', 'RT_end', 'RT_score', 'RT_length', 'tes_start', 'tes_end', 'gene_start', 'gene_end', 'gene_length', 'gene_rpkm', 'segment_start', 'segment_end', 'segment_score', 'segment_length', 'num_fused_seg', 'intersect_bp']
    reform_df[numeric_cols] = reform_df[numeric_cols].apply(pd.to_numeric)
    
    return reform_df


def filter_genes(df, assay):
    
    '''
    filters RT_segments to remove: 1.) segments that start with a score below 3,
    2.) genes that are shorter than 1kb, and 3.) genes with a gene_rpkm less than 0.25 
    or maybe less (assigned by assay)
    '''
    
    if assay == '0h' or assay == 'IR':
        rpkm = 0.25
    if assay == '2h':
        rpkm = 0.25
    if assay == '6h':
        rpkm = 0.25
    
    filtered_df = df[(df['segment_score'] > 2) & (df['gene_length'] >= 1000) & (df['gene_rpkm'] >= rpkm)]
    
    return filtered_df


def assign_rt_segments(df, assay):
    
    '''
    Parameters
    ----------
    df: takes reformatted RT_df output from reformat_df()
    
    Description
    -----------
    This function first filters and removes RT_segments based on desired segment
    score, gene_length, or gene_rpkm cutoffs [filter_genes()], it then reassigns
    duplicate RT_segments (those overlapping multiple genes) to the most upstream
    gene (**this can result in duplicates remaining in df if tes is shared between
    multiple unique genes)
    
    Variables
    ---------
    Cutoffs can be modified in filter_genes() function if desired
    '''
    #first take any gene_length or RPKM cutoffs and remove these rows
    #then assign the RT segment to the most upstream gene
    
    filt_df = filter_genes(df, assay)
    uniq_seg_counts = filt_df.segment_ID.value_counts()
    
    dup_seg = uniq_seg_counts[uniq_seg_counts > 1]
    dup_seg_list = list(dup_seg.index.values)
    
    dedup_df = filt_df.drop_duplicates(subset=['segment_ID'], keep=False, ignore_index=True)
    
    keep_df = pd.DataFrame()
    for seg_id in tqdm(dup_seg_list):
        eval_segs_df = filt_df[filt_df['segment_ID'].str.contains(seg_id)]
        if eval_segs_df['RT_strand'].iloc[0] == '-':
            tes_seg = eval_segs_df[eval_segs_df.RT_end == eval_segs_df.RT_end.max()]
        else:
            tes_seg = eval_segs_df[eval_segs_df.RT_start == eval_segs_df.RT_start.min()]
            
        keep_df = pd.concat([keep_df, tes_seg])

    rt_seg_df = pd.concat([dedup_df,keep_df], ignore_index=True)
    rt_seg_df_sorted = rt_seg_df.sort_values(by = ['RT_chr', 'RT_start'])
    
    return rt_seg_df_sorted


def main(bed_filepath, count_filepath, rpkm_filepath, sample, assay):
    
    '''
    Parameters
    ----------
    bed_filepath: [arg1] path to fused_dec bedfile (expects: .segments.fused_dec_simple.tes.intersect_test.bed)
    rpkm_filepath: [arg2] path to gene rpkm file (expects: make_counts_matrix_PEdata.sh output format)
    
    Description
    -----------
    Loads in fused_dec .bed output (intersected with tes) and sample genes file,
    filters and subsets both files, merges gene rpkm information with RT dataframe,
    refines RT_segments, and outputs refined RT_segment dataframe
    
    Variables
    ---------
    No variables to change in this main function.
    '''
    
    if bed_filepath is None:
        bed_filepath = "~/Desktop/shared_biosamples/read-through/fused_dec/HCT1160h4001a.hg38.gencode_29.spikeins_1.bin_250.segments.actual_RPKM.fused_dec.TESintersect.bed"
    
    if count_filepath is None:
        count_filepath = "~/Desktop/shared_biosamples/ALLcounts_ENCODE16CL_0h_genes.bed"
    
    if rpkm_filepath is None:
        rpkm_filepath = "~/Desktop/shared_biosamples/ALLrpkm_ENCODE16CL_0h_genes.bed"
    
    df = load_from_bed(bed_filepath)
    
    counts = pd.read_csv(count_filepath, sep='\t', engine='python')
    rpkm = pd.read_csv(rpkm_filepath, sep='\t', engine='python')
    
    rpkm_sub = subset_rpkm(counts, rpkm, sample)
    
    reformatted_df = reformat_df(df, rpkm_sub)
    
    refined_df = assign_rt_segments(reformatted_df, str(assay)) # sorted in this function
    
    rt_df = refined_df.rename(columns={"RT_chr":"#RT_chr"})
    
    
    bed_string = os.path.splitext(bed_filepath)[0]
    bed_output = bed_string.split(sep='.', maxsplit=1)[0]+".RTsegments.TESintersect.filtered.bed"
    delete_if_exists(bed_output)
    rt_df.to_csv(bed_output, sep='\t', header=True, index=False)


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_file',
                        help='Path to RT_seg file.')
    parser.add_argument('counts_file',
                        help='Path to counts matrix.')
    parser.add_argument('rpkm_file',
                        help='Path to RPKM matrix.')
    parser.add_argument('sample',
                        help='String denoting the sample name (i.e. A6730h4001a).')
    parser.add_argument('assay',
                        help='String denoting Bru assay. Required to appropriately assign rpkm cutoff (see filter_genes()).')
    
    args = parser.parse_args()

    main(args.bed_file, args.counts_file, args.rpkm_file, args.sample, args.assay)
    
    # bed_filepath = None
    # rpkm_filepath = None
    # main(bed_filepath, rpkm_filepath)
