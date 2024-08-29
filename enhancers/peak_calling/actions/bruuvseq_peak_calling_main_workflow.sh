#!/bin/bash


UV_GROUP="$1"; shift
SAMPLES=("$@")

STRANDS=("plus" "minus")

### activate conda environment

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV


### Run workflow

echo "-------------------------------"
echo
echo ${SAMPLES[@]}
echo $UV_GROUP
echo

echo "Separating .bam reads and calling peaks for individual samples..."
# separate reads for plus and minus strands and call peaks and summits for separate strands

CALL_SUMMITS=TRUE

# ***maybe put for loop in script and run that as array job

for SAMPLE in ${SAMPLES[@]}
do

	for STRAND in ${STRANDS[@]}
	do
	
	echo
	echo "Verifying loop: $SAMPLE $STRAND"
	echo
	
		. $ACTIONS_DIR/separate_strands.sh $SAMPLE $STRAND
		. $ACTIONS_DIR/macs2_callpeak.sh $SAMPLE $STRAND
	
	done
	
done


echo 
echo "-------------------------------"
echo 

echo "Calling peaks for replicate samples and defining reproducible peaks..."
# call peaks for replicate samples (no summits called) and get reproducible peaks between replicates

unset CALL_SUMMITS
CALL_SUMMITS=FALSE

for STRAND in ${STRANDS[@]}
do
	
	echo
	echo "Verifying loop: ${SAMPLES[@]} $STRAND"
	echo
	
	. $ACTIONS_DIR/macs2_callpeak_replicates.sh $STRAND "${SAMPLES[@]}"
	. $ACTIONS_DIR/reproducible_peaks.sh $STRAND "${SAMPLES[@]}"
	
done


echo 
echo "-------------------------------"
echo 

echo "Recombining strands..."
# combine plus and minus strand results

. $ACTIONS_DIR/combine_strands.sh "${SAMPLES[@]}"


echo 
echo "-------------------------------"
echo 

echo "Calculating FRIP score..."
# calculate frip

. $ACTIONS_DIR/frip.sh "${SAMPLES[@]}"


echo 
echo "-------------------------------"
echo 

echo "BruUV-seq peak calling complete!"