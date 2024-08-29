#!/bin/bash

### description


### activate conda environment

CONDA_ENV=/nfs/turbo/radonc-ljungman-turbo/envs/p39							# path to conda environment

conda activate $CONDA_ENV


### General variables

WORKING_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/bipeaks/working
SAMPLE_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/bipeaks/sample		# location of reproducible peak files and cCRE files
ACTION_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/bipeaks/actions

# SAMPLES=("A673" "Caco2" "Calu3" "GM12878" "HCT116" "HepG2" "HMEC" "HUVEC"
		 # "IMR90" "K562" "MCF10A" "MCF7" "OCILY7" "panc1" "PC3" "PC9"
		 # ) # technically current peak_calling workflow outputs files as "A673UV.reproducible.peak" so would need to modify this to use in tandem

SAMPLES=("A673"
		 "HCT116"
		 "PC3")

DIST=600		# bp distance allowed between bidirectional enhancer peak pairs


### Run workflow

echo "-------------------------------"
echo 

PEAK_DIR=$WORKING_DIR/peak_dir_output

mkdir -p $PEAK_DIR

echo "Pairing bidirectional peaks..."
# define bipeak pairs from reproducible peaks (expects $SAMPLE.reproducible.peaks.bed filename convention)

# ***expects input file list in dir; BASEDIR and OUTDIR specification is optional -- maybe modify this

. $ACTION_DIR/bipeaks_AM.sh $SAMPLE_DIR $WORKING_DIR/peak_dir_output ${SAMPLES[@]}


echo 
echo "-------------------------------"
echo 


echo "Extract bidirectional peak pairs that meet distance threshold..."
# identify bipeaks and unipeaks according to threshold (600bp) and creates bipeak, unipeak, and allpeak files

for SAMPLE in ${SAMPLES[@]}
do

	echo $SAMPLE
	
	# python $ACTION_DIR/generate_peak_files_v2.py $PEAK_DIR/$SAMPLE.bipeaks.bed --dist $DIST		# v2 is really slow (probably why I used original despite warnings for deprecated function)
	python $ACTION_DIR/generate_peak_files.py $PEAK_DIR/$SAMPLE.bipeaks.bed --dist $DIST
	
done


echo 
echo "-------------------------------"
echo 


echo "Obtaining peaks overlapping enhancers..."
# filter reproducible peaks for eRNA peaks

ENH_DIR=$WORKING_DIR/enhancer_peaks

mkdir -p $ENH_DIR

for SAMPLE in ${SAMPLES[@]}
do
	
	echo $SAMPLE
	
	# intersect peaks with dELS cCREs
	bedtools intersect -wo -a $SAMPLE_DIR/$SAMPLE.reproducible.peaks.bed -b $SAMPLE_DIR/$SAMPLE-dELS.bed | \
		
	# remove peaks overlapping pELS cCREs
	bedtools intersect -v -a - -b $SAMPLE_DIR/$SAMPLE-pELS.bed | \
	
	# collapse output to remove duplicate entries
	bedtools groupby -i - -g 1-10 -c 14 -o count > $ENH_DIR/$SAMPLE.enhancer_peaks.bed
	
done


echo 
echo "-------------------------------"
echo 


echo "Filtering bidirectional peaks for those overlapping enhancers..."
# filter files for eRNA peaks (those that overlap dELS and do not overlap pELS cCREs).

. $ACTION_DIR/bipeak_dels_filter.sh ${SAMPLES[@]}
. $ACTION_DIR/unipeak_dels_filter.sh ${SAMPLES[@]}
. $ACTION_DIR/allpeak_dels_filter.sh ${SAMPLES[@]}


echo 
echo "-------------------------------"
echo 


echo "Finding consensus peak regions across cell lines..."
# find consensus peaks across all samples (overlapping by 90%) and merged into a single dataframe

OL_DIR=$WORKING_DIR/overlap_files

mkdir -p $OL_DIR

for SAMPLE in ${SAMPLES[@]}
do

	COMMAND="bedtools intersect -wao -f 0.90 -F 0.90 -e"
	echo $SAMPLE
	COMMAND+=" -a $PEAK_DIR/$SAMPLE.unique_peaks_$DIST.dELS_intersect.bed -b"
	
	for S in ${SAMPLES[@]}
	do
		
		if [ $S == $SAMPLE ]; then
			continue
		
		else
			echo $S
			COMMAND+=" $PEAK_DIR/$S.unique_peaks_$DIST.dELS_intersect.bed"
		fi
		
	done

	COMMAND+=" > ${OL_DIR}/${SAMPLE}-all.peak_overlaps.bed"
	echo $COMMAND
	eval $COMMAND
	
done

# this script expects input file naming convention established above (${SAMPLE}-all.peak_overlaps.bed); matrix output defaults to include all DPCL cell lines even if peak files are not supplied
python $ACTION_DIR/overlapping_peak_files_v2.py $OL_DIR

mv BruUV_eRNA_16_cell_lines.tsv $OL_DIR
mv BruUV_eRNA_16_cell_lines.bed $OL_DIR

echo 
echo "-------------------------------"
echo 


### will add this in once revisions are made
# echo "Classifying enhancer peaks..."

# #get counts for peak regions above

# # classify peaks
# python $ACTION_DIR/classify_peaks.py



echo 
echo "-------------------------------"
echo 

echo "Enhancer analysis complete!"





