#!/bin/bash

# This script should only be run once per sample. This will copy segments.actual_RPKM files from the develop samples folder to the working directory and 
# unzips the files for use in downstream steps.

#q    require $SAMPLE $DEV_SAMPLES_DIR $GB_RAM
#q    require $BIN_SIZE $ANNOTATION $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $DEV_SAMPLES_DIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $DEV_SAMPLES_DIR
#PBS  -l  mem=$GB_RAM\gb

# cp sample segments data (actual_RPKM only)
echo "copying file: $SAMPLE_DIR/$SAMPLE/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.bed.bgz"
cp $DEV_SAMPLES_DIR/$SAMPLE/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.bed.bgz $OUTDIR

# unzip segments file
echo "unzipping files"
gunzip -dc $OUTDIR/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.bed.bgz > $OUTDIR/$SAMPLE.$ANNOTATION.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.bed
