#!/bin/bash

# This script run bedtools intersect to find fused_dec segments that overlap a TES. A_file = TES_annotation so that each row in the output cooresponds to a gene.

#q    require $SAMPLE $GB_RAM $TES_ANNOT
#q    require $BIN_SIZE $ANNOTATION $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb

echo "intersect segments with TESs"
bedtools intersect -s -wao -a $TES_ANNOT -b $OUTDIR/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.fused_dec.bed > $OUTDIR/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.fused_dec.TESintersect.bed

