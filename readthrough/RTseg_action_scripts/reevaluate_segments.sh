#!/bin/bash

# This script uzips count_bed outputs of eval regions and runs the python script that evaluates these segments and reassigns RT_seg coordinates.

#q    require $SAMPLE $ACTIONS_DIR $GB_RAM
#q    require $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb


# unzip count_bed output
gunzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.startsRPKM.bed.bgz > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.startsRPKM.bed
gunzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.endsRPKM.bed.bgz > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.endsRPKM.bed
gunzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genesRPKM.bed.bgz > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genesRPKM.bed

# reevaluate ambiguous segment and add back to full df
conda run -n bruseq5 python $ACTIONS_DIR/readthru/evaluate_ambiguous_segs_v2.py \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.startsRPKM.bed \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.endsRPKM.bed \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genesRPKM.bed \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.noeval.bed \
$OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.eval.bed
