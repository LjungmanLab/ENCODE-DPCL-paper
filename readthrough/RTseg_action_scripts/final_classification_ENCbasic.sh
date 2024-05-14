#!/bin/bash

# This script unzips the count_bed output from the previous step and runs the python script that performs the final classification of all RT_segments
# and formats the final output file.

#q    require $SAMPLE $ACTIONS_DIR $GB_RAM
#q    require $OUTDIR $ANNOTATION $RPKM_6h $RPKM_0h

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb


# re-filter and classify updated segments
conda run -n bruseq5 python $ACTIONS_DIR/readthru/final_filter_and_classification_v2_ENCbasic.py \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed \
$SAMPLE \
$RPKM_0h \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.closestDS.bed \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.reverseOL.bed \
--exon_expression_file $RPKM_6h \