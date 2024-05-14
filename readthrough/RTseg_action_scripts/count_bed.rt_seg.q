
###########

$GENOME        hg38
$ANNOTATION    hg38.gencode_29
$CUSTOM_GENOME spikeins_1
$STRANDED      True

invoke $ACTIONS_DIR/file_info/genome_info.q


invoke $ACTIONS_DIR/file_info/genome_info.q
invoke $ACTIONS_DIR/file_info/sample_ann_info.q
invoke $ACTIONS_DIR/file_info/sample_ann_info.q

$COV_TYPE generic
$FEAT_BED_IN  $IN_BED
$FEAT_BED_OUT $OUT_BED
$GB_RAM 2

$S_OPTION  $STRANDED ? "-s" : ""

qsub $ACTIONS_DIR/count/feature_coverage.sh










# recreation of count_bed.q for use in RT_segment workflow
# IN_BED, OUT_BED, and GB_RAM are specified in main q script

# $STRANDED		True

# invoke			$ACTIONS_DIR/file_info/genome_info.q


# invoke $ACTIONS_DIR/file_info/genome_info.q
# invoke $ACTIONS_DIR/file_info/sample_ann_info.q
# invoke $ACTIONS_DIR/file_info/sample_ann_info.q

# $COV_TYPE		generic 
# $FEAT_BED_IN	$IN_BED 
# $FEAT_BED_OUT	$OUT_BED

# $S_OPTION  $STRANDED ? "-s" : ""

# qsub $ACTIONS_DIR/count/feature_coverage.sh
