#!/bin/bash

FILEDIR="$1"; shift
OUTDIR="$1"; shift
SAMPLES=("$@")

if [[ "$FILEDIR" == "" ]]
then
	echo "FILEDIR not specified, defaulting to: ./peak_pairs"
	FILEDIR="/peak pairs"
fi

if [[ "$OUTDIR" == "" ]]
then
	echo "OUTDIR not specified, defaulting to: ./bidir_peaks_output"
	OUTDIR="/bidir_peaks_output"
fi

mkdir -p $OUTDIR

# function to find closest peak on opposite strand
function bipeaks_5p {

	INFILE="$1"
	OUTPRE="$2"
#	MAX_DIST=1000
# | sort -k1,1 -k2,2n)


	# note we are adding 2 because we will be adding 2 cols before bedtools
	A_COLNUM=$(gawk 'BEGIN{FS=OFS="\t"}(FNR==1){print NF+2}' $INFILE)
	B_COLNUM=$A_COLNUM


	bedtools closest -S -D a -nonamecheck \
		-a <( gawk \
			'BEGIN{
				FS=OFS="\t"
			}($1~/^chr[^_]+$/){
				start=$2
				end=$3
				if($6=="+"){
					$2=$2
					$3=$2+1
				}else{
					$2=$3-1
					$3=$3
				}
				print $0,start,end,"@@@@"
			}' $INFILE \
			| sort -k 1,1 -k 2,2n ) \
		-b <( gawk \
			'BEGIN{
				FS=OFS="\t"
			}($1~/^chr[^_]*$/){
				start=$2
				end=$3
				if($6=="+"){
					$2=$2
					$3=$2+1
				}else{
					$2=$3-1
					$3=$3
				}
				print $0,start,end,"@@@@"
			}' $INFILE \
			| sort -k 1,1 -k 2,2n ) \
		| gawk 'BEGIN{FS=OFS="\t"}{printf $1; for(i=2;i<=NF;i++){ if($i != "@@@@"){printf OFS $i}} printf ORS}' \
		| gawk \
		-v a_colnum=$A_COLNUM \
		-v b_colnum=$B_COLNUM \
		-v cell_line=$CL \
		'# min value, including tie
		function find_min(num1, num2){
		   if(num1 <= num2){
				return num1
			}else{
				return num2
			}
		}
		# max value, including tie
		function find_max(num1, num2){
			if(num1 >= num2){
				return num1
			}else{
				return num2
			}
		}
		BEGIN{
			FS=OFS="\t"
			distance_idx = a_colnum + b_colnum + 1
			ct=0
		}{
			ct++
			d = $distance_idx < 0 ? -($distance_idx) : $distance_idx

			# define bipeak coords using min and max coords of the 2 peaks
			min_coord = find_min( \
				find_min( \
					find_min( \
						$(a_colnum - 1),$a_colnum \
					),$(a_colnum + b_colnum - 1) \
				),$(a_colnum + b_colnum) \
			)
			max_coord = find_max( \
				find_max( \
					find_max( \
						$(a_colnum - 1),$a_colnum \
					),$(a_colnum + b_colnum - 1) \
				),$(a_colnum + b_colnum) \
			)


			# main bipeak BED definition
			outline = $1 OFS min_coord OFS max_coord OFS cell_line "UV.bipeak_" ct OFS 0 OFS "."
			# add the midpoint of 5p coords
			outline = outline OFS int($2+(($(a_colnum+2)-$2)/2))

			A_peak_cols = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6
			if(a_colnum > 6 ){
				for(i=7;i<=a_colnum;i++){
					A_peak_cols = A_peak_cols OFS $i
				}
			}
			B_peak_cols = $(a_colnum + 1) OFS $(a_colnum + 2) OFS $(a_colnum + 3) OFS $(a_colnum + 4) OFS $(a_colnum + 5) OFS $(a_colnum + 6)
			if(b_colnum > 6 ){
				for(i=7;i<=b_colnum;i++){
					B_peak_cols = B_peak_cols OFS $(a_colnum + i)
				}
			}
			# make sure plus strand peak is always the A peak
			if($6=="+"){
				outline = outline OFS A_peak_cols
				outline = outline OFS B_peak_cols
			}else{
				outline = outline OFS B_peak_cols
				outline = outline OFS A_peak_cols
			}

			if(max_coord == $a_colnum && min_coord==$(a_colnum+b_colnum-1)){
				if($6=="+"){
					class = "divergent"
				}else{
					class = "convergent"
				}
			}else if(max_coord == $a_colnum && min_coord==$(a_colnum-1)){
				if($6=="+"){
					class = "AcoversB"
				}else{
					class = "BcoversA"
				}
			}else if(max_coord == $(a_colnum+b_colnum) && min_coord==$(a_colnum+b_colnum-1)){
				if($6=="+"){
					class = "BcoversA"
				}else{
					class = "AcoversB"
				}
			}else if(max_coord == $(a_colnum+b_colnum) && min_coord==$(a_colnum-1)){
				if($6=="+"){
					class = "convergent"
				}else{
					class = "divergent"
				}
			}else{
				class = "error"
			}
			outline = outline OFS d OFS $distance_idx OFS class
			print outline
		}' \
	| sort -k1,1 -k2,2n \
	> $OUTPRE.bipeaks.bed

	# cat $OUTFILE_PREFIX.txt \
	# | gawk 'BEGIN{FS=OFS="\t"}{
		# print $1,$7,$7+1,$4,$5,$6,$7,$7+1,"0,0,0"
		# }' \
		# | sort -k1,1 -k2,2n \
	# > $OUTFILE_PREFIX.bed9 

	# bedToBigBed $OUTFILE_PREFIX.bed9 hg38.chromSizes $OUTFILE_PREFIX.bb

}

for CL in ${SAMPLES[@]}
	do echo $CL
	IN_PEAKS_FILE=$FILEDIR/$CL.reproducible.peaks.bed
	if [[ -f "$IN_PEAKS_FILE" ]]
	then
		bipeaks_5p $IN_PEAKS_FILE $OUTDIR/$CL
	else
		echo "Peak file $IN_PEAKS_FILE does not exist. Skipping..."
	fi

done

