#!/bin/bash

# This script unzips the previous count_bed output, and runs the python script that determines the segments that can be reevaluated and
# outputs the regions needed for evaluation (1kb into start of RT_seg, 1kb before downstream TSS, 1kb after downstream TSS). These outputs
# are then zipped for count_bed.

#q    require $SAMPLE $ACTIONS_DIR $GB_RAM
#q    require $RPKM_6h $ANNOTATION $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb


# unzip count_bed output
gunzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.bed.bgz > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.bed

# assign preliminary RT_classes and get eval regions of ambiguous segments
conda run -n bruseq5 python $ACTIONS_DIR/readthru/get_rt_eval_regions_v3_ENCbasic.py $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.bed $SAMPLE --exon_expression_file $RPKM_6h

# zip python output for count_bed
gzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.starts.bed > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.starts.bed.bgz
gzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.ends.bed > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.ends.bed.bgz
gzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genes.bed > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.genes.bed.bgz
