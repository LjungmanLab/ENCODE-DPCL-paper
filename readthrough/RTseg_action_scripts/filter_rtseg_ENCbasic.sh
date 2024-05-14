#!/bin/bash

# This script runs the python script that performs filtering and refining of readthrough segments. Inputs include a .bed file containing RT_segments (fused_dec
# segments that intersect a TES) a file that contains the RPKM information for the sample () and the sample name.
# The output is a .filtered.bed file that follows the naming of the input file (defined in python main func). More information can be found within python script.

#q    require $SAMPLE $ACTIONS_DIR $GB_RAM $COUNTS_0h $BRU_TP
#q    require $RPKM_0h $BIN_SIZE $ANNOTATION $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb


# filter and refine rt segments
echo "filter and refine RT_segments"
conda run -n bruseq5 python $ACTIONS_DIR/readthru/filter_and_refine_rt_segments_v2_ENCbasic.py \
$OUTDIR/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.fused_dec.TESintersect.bed \
$COUNTS_0h \
$RPKM_0h \
$SAMPLE \
$BRU_TP
