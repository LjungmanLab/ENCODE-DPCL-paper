# readthru script description

# General variables

$PROJECT				ENCODE
$SUB_PROJECT				16_cellline_0h_run4
$OUTDIR            			$MYPROJECTS_DIR/readthru/$PROJECT/$SUB_PROJECT
$GB_RAM					2

$GENOME					hg38
$ANNOTATION				hg38.gencode_29
$CUSTOM_GENOME				spikeins_1
$BRU_TP					0h
$DEV_SAMPLES_DIR			/treehouse/wilsonte_lab/radonc_ljungman_turbo/Ljungman5-develop/samples

$BIN_SIZE				250
$VERSION				v8
$GENE_ANNOT				/home/wilsonte_lab/clubhouse/usr/mcshanea/ENCODE/readthru/fused_dec/gencode.v29.basic.annotation_genes_final_sorted.bed
$TES_ANNOT				/home/wilsonte_lab/clubhouse/usr/mcshanea/ENCODE/readthru/fused_dec/gencode.v29.basic.annotation_tes_final_sorted.bed
$COUNTS_0h				/home/wilsonte_lab/clubhouse/usr/mcshanea/ENCODE/readthru/fused_dec/ENCbasic_0h_gene_counts/ALLcounts_ENCODE16CL_0h_genes.bed
$RPKM_0h				/home/wilsonte_lab/clubhouse/usr/mcshanea/ENCODE/readthru/fused_dec/ENCbasic_0h_gene_counts/ALLrpkm_ENCODE16CL_0h_genes.bed
$RPKM_6h				/home/wilsonte_lab/clubhouse/usr/mcshanea/ENCODE/readthru/fused_dec/ENCbasic_6h_exon_counts/ALLrpkm_ENCODE16CL_6h_exons.bed

$SAMPLES				HUVEC0h4002a

# $SAMPLES				A6730h4001a A6730h4002a Caco20h4001a Caco20h4002a Calu30h4001a Calu30h4002a GM128780h4001a GM128780h4002a
# $SAMPLES $SAMPLES			HCT1160h4001a HCT1160h4002a HepG20h4005a HepG20h4006a HMEC0h4001a HMEC0h4002a HUVEC0h4001a HUVEC0h4002a
# $SAMPLES $SAMPLES			IMR900h4003a IMR900h4004a K5620h4004a K5620h4005a MCF10A0h4001a MCF10A0h4002a MCF70h4001a MCF70h4002a
# $SAMPLES $SAMPLES			OCILY70h4001a OCILY70h4002a panc10h4001a panc10h4002a PC30h4001a PC30h4002a PC90h4001a PC90h4002a

invoke get_rt_segs.q $SAMPLE $SAMPLES

<file name="get_rt_segs.q">


# copy files
$STEP_NAME				copy

thread cp
qsub $ACTIONS_DIR/readthru/copy_seg_files.sh


# run fused_dec python script
$STEP_NAME				fuse_seg

thread fused_dec cp
qsub $ACTIONS_DIR/readthru/create_fused_dec_seg.sh


# run bedtools intersect of RT_segs and gene TESs
$STEP_NAME				tes_int

thread intersect fused_dec
qsub $ACTIONS_DIR/readthru/segment_tes_intersect.sh


# run filtering and refining python script
$STEP_NAME				filter

thread filter intersect
qsub $ACTIONS_DIR/readthru/filter_rtseg_ENCbasic.sh


# run bedtools intersect of RT_segs and genes
$STEP_NAME				ct_ol_1
thread overlap1 filter

$IN_FILE	$OUTDIR/$SAMPLE.RTsegments.TESintersect.filtered.bed
$OUT_FILE	$OUTDIR/$SAMPLE.RTsegments.TESintersect.geneOL.bed

qsub $ACTIONS_DIR/readthru/count_overlapping_genes.sh


# run count_bed to get RT_seg counts
thread count_bed1 overlap1

$IN_BED		$OUTDIR/$SAMPLE.RTsegments.TESintersect.geneOL.bed.bgz
$OUT_BED	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.bed.bgz

invoke $ACTIONS_DIR/readthru/count_bed.rt_seg.q


# run python script to get eval regions of ambiguous RT_segs
$STEP_NAME				eval_reg

thread eval_reg count_bed1
qsub $ACTIONS_DIR/readthru/get_eval_regions_ENCbasic.sh


# run count_bed to get counts of eval regions
thread count_bed2 eval_reg

$IN_BED		$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.starts.bed.bgz
$OUT_BED	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.startsRPKM.bed.bgz

invoke $ACTIONS_DIR/readthru/count_bed.rt_seg.q

$IN_BED		$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.ends.bed.bgz
$OUT_BED	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.endsRPKM.bed.bgz

invoke $ACTIONS_DIR/readthru/count_bed.rt_seg.q

$IN_BED		$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genes.bed.bgz
$OUT_BED	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genesRPKM.bed.bgz

invoke $ACTIONS_DIR/readthru/count_bed.rt_seg.q


# run python script to reevaluate ambiguous RT_segs
$STEP_NAME				reeval_seg

thread evaluate count_bed2
qsub $ACTIONS_DIR/readthru/reevaluate_segments.sh


# re-run bedtools intersect for final RT_segs and genes
$STEP_NAME				ct_ol_2

thread overlap2 evaluate

$IN_FILE	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixed.bed
$OUT_FILE	$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedOL.bed

qsub $ACTIONS_DIR/readthru/count_overlapping_genes.sh


# re-run count_bed to get counts of final RT_segs
thread count_bed3 overlap2

$IN_BED      $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedOL.bed.bgz
$OUT_BED     $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed.bgz
 
invoke $ACTIONS_DIR/readthru/count_bed.rt_seg.q

# get the closest downstream genes and overlapping reverse strand genes for final classification
$STEP_NAME				close_and_rev

thread close_rev count_bed3
qsub $ACTIONS_DIR/readthru/rev_ol_and_closest_ds.sh


# run final processing python script
$STEP_NAME				final

thread final close_rev
qsub $ACTIONS_DIR/readthru/final_classification_ENCbasic.sh


</file>



