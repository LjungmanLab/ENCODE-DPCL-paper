#!/bin/bash

SAMPLE_ARR=("$@")

for SAMPLE in ${SAMPLE_ARR[@]}
do

	FILE_NAME=$PEAK_DIR/$SAMPLE.unique_unipeaks_$DIST.bed

	echo $FILE_NAME

	Rscript $ACTION_DIR/unipeak_dels_filter.R $FILE_NAME $ENH_DIR/$SAMPLE.enhancer_peaks.bed

done
