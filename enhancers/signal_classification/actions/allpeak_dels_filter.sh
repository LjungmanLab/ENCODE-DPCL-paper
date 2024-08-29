#!/bin/bash

SAMPLE_ARR=("$@")

for SAMPLE in ${SAMPLE_ARR[@]}
do

	FILE_NAME=$PEAK_DIR/$SAMPLE.unique_peaks_$DIST.bed
	BIPEAKS=$PEAK_DIR/$SAMPLE.unique_bipeaks_$DIST.dELS_intersect.bed
	UNIPEAKS=$PEAK_DIR/$SAMPLE.unique_unipeaks_$DIST.dELS_intersect.bed

	echo $FILE_NAME

	Rscript $ACTION_DIR/allpeak_dels_filter.R $FILE_NAME $BIPEAKS $UNIPEAKS

done
