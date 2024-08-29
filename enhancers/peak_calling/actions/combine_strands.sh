#!/bin/bash

### description

# script written by Brian Magnuson; script modified by Ariel McShane (20240819)


SAMPLES=("$@")

# CHROMS_LIST=$(echo "$CHROMS_STD" | tr '\n' ',' | sed 's/ /,/g')
CHROMS_LIST="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

function combine_strands {
	
	OUT_PRE=$1
	OUT_EXT=$2
	PLUS_FILE=$3
	MINUS_FILE=$4

	# combine strands AND filter out non-canonical and spike-in chroms
	cat \
		<(gawk 'BEGIN{FS=OFS="\t"}{$6="+"; print $0}' $PLUS_FILE) \
		<(gawk 'BEGIN{FS=OFS="\t"}{$6="-"; print $0}' $MINUS_FILE) \
        | gawk \
                -v canonical_chroms_string="$CHROMS_LIST" \
                'BEGIN{
                        FS=OFS="\t"
                        split(canonical_chroms_string,cc,",")
                        chr_regex = "^" cc[1] "$"
                        for(i=2;i<=length(cc);i++){
                                chr_regex = chr_regex "|^" cc[i] "$"
                        }
                }($1~chr_regex){
			score = $5 > 1000 ? 1000 : $5  # max score value for UCSC standard is 1000, so saturate this value
			$5 = score
			print $0
                }' \
	| sort -k1,1 -k2,2n \
	> $OUT_PRE.$OUT_EXT

	# gzip
	gzip -c $OUT_PRE.$OUT_EXT > $OUT_PRE.$OUT_EXT.gz
	cp $OUT_PRE.$OUT_EXT.gz $OUT_PRE.ENCODE.$OUT_EXT.gz

	
	#make bed9 and bigbed
	cat $OUT_PRE.$OUT_EXT \
	| gawk \
		'BEGIN{
			FS=OFS="\t"
		}{
			color = $6=="+" ? "30,144,255" : "255,30,143"
			score = $5 > 1000 ? 1000 : $5  # max score value for UCSC standard is 1000, so saturate this value
			print $1,$2,$3,$4,score,$6,$2,$3,color
		}' \
	> $OUT_PRE.bed9

	bedToBigBed $OUT_PRE.bed9 $CHROM_SIZES $OUT_PRE.bb

	# ENCODE related
	gzip -c $OUT_PRE.bed9 > $OUT_PRE.ENCODE.bed.gz  # ENCODE wants a gzipped bed file
	cp $OUT_PRE.bb $OUT_PRE.ENCODE.bb # bb are identical, but we want ENCODE in the name for consistency
}

# ### load conda env

# if [ -f "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh" ]
# then
	# . "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh"
# fi

# OLD_ENV=""
# if [[ "$CONDA_DEFAULT_ENV" != "bruseq5" ]]
# then
	# if [[ "$CONDA_DEFAULT_ENV" == "base" ]]
	# then
		# conda deactivate
	# elif [[ "$CONDA_DEFAULT_ENV" != "" ]]
	# then
		# OLD_ENV=$CONDA_DEFAULT_ENV
		# conda deactivate
	# fi
	# conda activate bruseq5
# fi

#######################################
# combined strands and make final files

# unfiltered pooled replicates peaks
POOLEDREPS_PEAKS_PLUS=$PROJECT_DIR/$UV_GROUP.plus/$UV_GROUP.plus_peaks.narrowPeak
POOLEDREPS_PEAKS_MINUS=$PROJECT_DIR/$UV_GROUP.minus/$UV_GROUP.minus_peaks.narrowPeak
combine_strands $PROJECT_DIR/$UV_GROUP.pooled_replicates.peaks narrowPeak $POOLEDREPS_PEAKS_PLUS $POOLEDREPS_PEAKS_MINUS

# reproducible peaks (filtered pooled replicates peaks)
# note: output ext is narrowPeak even though input are bed files -- bed are narrowPeak compatible
REPRODUCIBLE_PEAKS_PLUS=$PROJECT_DIR/$UV_GROUP.plus.reproducible.peaks.bed
REPRODUCIBLE_PEAKS_MINUS=$PROJECT_DIR/$UV_GROUP.minus.reproducible.peaks.bed
combine_strands $PROJECT_DIR/$UV_GROUP.reproducible.peaks narrowPeak $REPRODUCIBLE_PEAKS_PLUS $REPRODUCIBLE_PEAKS_MINUS

# peaks for indiv samples
for SAMPLE in ${SAMPLES[@]}
do
	NARROWPEAK_PLUS_FILE=$PROJECT_DIR/$SAMPLE.plus.summits/$SAMPLE.plus_peaks.narrowPeak
	NARROWPEAK_MINUS_FILE=$PROJECT_DIR/$SAMPLE.minus.summits/$SAMPLE.minus_peaks.narrowPeak
	combine_strands $PROJECT_DIR/$SAMPLE.peaks narrowPeak $NARROWPEAK_PLUS_FILE $NARROWPEAK_MINUS_FILE 
done
#######################################

# ### unload conda env
# conda deactivate
# # reactivate old env (if applicable)
# if [[ "$OLD_ENV" != "" ]]
# then
	# conda activate $OLD_ENV
# fi



