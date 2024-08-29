#!/bin/bash

### description

# script written by Brian Magnuson; script modified by Ariel McShane (20240819)


STRAND="$1"; shift
SAMPLES=("$@")

REPS_PEAKS=$PROJECT_DIR/$UV_GROUP.$STRAND/$UV_GROUP.${STRAND}_peaks.narrowPeak
unset SUMMIT_FILES
SUMMIT_FILES=($(for SAMPLE in ${SAMPLES[@]}; do echo -n " $PROJECT_DIR/$SAMPLE.$STRAND.summits/$SAMPLE.${STRAND}_summits.bed";done))
OUTFILE=$PROJECT_DIR/$UV_GROUP.$STRAND.reproducible.peaks.bed


function reproducible_peaks {

	REPA_SUMMITS_TMP=${SUMMIT_FILES[0]}.tmp
	REPB_SUMMITS_TMP=${SUMMIT_FILES[1]}.tmp
	
	lasti=$(echo "${#SUMMIT_FILES[@]} - 1" | bc -l)
	echo "${SUMMIT_FILES[0]}"

	for i in $(seq 1 $lasti)
	do
		echo "${SUMMIT_FILES[$i]}"
		if [[ $i == 1 ]]
		then
			A_TMP=${SUMMIT_FILES[0]}.tmp
			B_TMP=${SUMMIT_FILES[1]}.tmp
			
			bedtools intersect \
				-nonamecheck -u \
				-a $REPS_PEAKS  \
				-b <(gawk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,0,".",$5}' ${SUMMIT_FILES[0]})  \
			> $A_TMP
			
			bedtools intersect \
				-nonamecheck -u \
				-a $REPS_PEAKS  \
				-b <(gawk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,0,".",$5}' ${SUMMIT_FILES[1]})  \
			> $B_TMP

		else
			lasti=$(echo "$i - 1" | bc -l)
			A_TMP=$OUTFILE.tmp.$lasti
			B_TMP=${SUMMIT_FILES[$i]}.tmp		

			bedtools intersect \
				-nonamecheck -u \
				-a $REPS_PEAKS  \
				-b <(gawk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,0,".",$5}' ${SUMMIT_FILES[$i]})  \
			> $B_TMP

		fi
		
		bedtools intersect \
			-nonamecheck -u \
			-a $A_TMP \
			-b $B_TMP \
		> $OUTFILE.tmp.$i

	done

	
	cp $OUTFILE.tmp.$lasti $OUTFILE

}

#######################################
# filter reproducible peaks
reproducible_peaks
#######################################
