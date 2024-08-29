#!/bin/bash


unset SAMPLES_ARR
SAMPLES_ARR=("$@")

REPRODUCIBLE_PEAKS_PLUS=$PROJECT_DIR/$UV_GROUP.plus.reproducible.peaks.bed
REPRODUCIBLE_PEAKS_MINUS=$PROJECT_DIR/$UV_GROUP.minus.reproducible.peaks.bed
OUTFILE=$PROJECT_DIR/$UV_GROUP.frip.txt


function frip {

	# create frip report
	echo -e "sample_pair_name\tsample\treplicate\tstrand\tkey\tvalue" > $OUTFILE
	
	for i in ${!SAMPLES_ARR[@]}
	do
		SAMPLE=${SAMPLES_ARR[$i]}
		PLUS_READS=$(samtools view -@ 16 -c $PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.plus.bam)
		MINUS_READS=$(samtools view -@ 16 -c $PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.minus.bam)
		PLUS_PEAKS_READS=$(samtools view -@ 16 -c -L $REPRODUCIBLE_PEAKS_PLUS $PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.plus.bam)
		MINUS_PEAKS_READS=$(samtools view -@ 16 -c -L $REPRODUCIBLE_PEAKS_MINUS $PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.minus.bam)

		FRIP_P=$(gawk 'BEGIN{print '$PLUS_PEAKS_READS'/'$PLUS_READS'}')
		FRIP_M=$(gawk 'BEGIN{print '$MINUS_PEAKS_READS'/'$MINUS_READS'}')

		echo -e "$UV_GROUP\t$SAMPLE\t1\tplus\ttotal_reads\t$PLUS_READS" >> $OUTFILE
		echo -e "$UV_GROUP\t$SAMPLE\t1\tplus\treads_in_peaks\t$PLUS_PEAKS_READS" >> $OUTFILE
		echo -e "$UV_GROUP\t$SAMPLE\t1\tplus\tfrip\t$FRIP_P" >> $OUTFILE

		echo -e "$UV_GROUP\t$SAMPLE\t1\tminus\ttotal_reads\t$MINUS_READS" >> $OUTFILE
		echo -e "$UV_GROUP\t$SAMPLE\t1\tminus\treads_in_peaks\t$MINUS_PEAKS_READS" >> $OUTFILE
		echo -e "$UV_GROUP\t$SAMPLE\t1\tminus\tfrip\t$FRIP_M" >> $OUTFILE
	done
}

#######################################
# combined strands and make final files
frip
#######################################



