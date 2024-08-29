#!/bin/bash

# recreation of pipeline `feature_coverage.sh` script from pipeline for slurm compatibility -- AM 20240719
# script calculates RPKM for sorted (b)gzipped input BED features in FEAT_BED_IN


### argument variables
SAMPLE=$1
FEAT_BED_IN=$2
FEAT_BED_OUT=$3


# ### load modules
# echo "Loading modules..."

# ml load Bioinformatics
# ml load bedtools2/2.30.0-svcfwbm
# ml load Rtidyverse/4.4.0

### activate conda environment

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV


### report program versions
echo
echo "Program versions:"
which bedtools
which bgzip
which Rscript
echo


echo "Calculating coverage for $SAMPLE"
echo


### set additional variables

## bedtools command options
S_OPTION=$([[ "$STRANDED" == True ]] && echo "-s" || echo "")

## sample directory and prefix
SAMPLE_DIR=/nfs/turbo/radonc-ljungman-turbo/Ljungman5/samples/$SAMPLE
SAMPLE_ROOT=$SAMPLE_DIR/$SAMPLE

# SMP_GEN_ROOT=$([[ "$CUSTOM_GENOME" == "" ]] && echo "$SAMPLE_ROOT.$ANNOTATION" || echo "$SAMPLE_ROOT.$ANNOTATION.$CUSTOM_GENOME")

## log directory and prefix
LOG_DIR=$SAMPLE_DIR/logs
LOG_PREFIX=$LOG_DIR/$SAMPLE

# STATS_FILE=$LOG_PREFIX.stats.txt
# STATS_FILE=$([[ "$CUSTOM_GENOME" == "" ]] && echo "$LOG_PREFIX.$ANNOTATION.stats.txt" || echo "$LOG_PREFIX.$ANNOTATION.$CUSTOM_GENOME.stats.txt")
STATS_FILE=$WORKING_DIR/ENCODE16CL_nMapUniq_values.txt
# STATS_FILE=$LOG_PREFIX.hg38.gencode_29.spikeins_1.stats.txt

## genome base coverage file
# BASE_COV_FILE=$SMP_GEN_ROOT.baseCoverage.bgz
BASE_COV_FILE=$SAMPLE_ROOT.hg38.gencode_29.spikeins_1.baseCoverage.bgz

## set script path
PFX=$ACTIONS_DIR
# PFX=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/actions

## get mapped read count for RPKM calculation
# nMapUniq=`awk '/^nMapUniq=/{ print substr($0,10)}' $STATS_FILE`
nMapUniq=$(awk -v s="$SAMPLE" '$1==s {print $2}' $STATS_FILE)

### double check correct variables were assigned

echo "Check assigned variables:"
echo "Sample: $SAMPLE"
echo "Stats file: $STATS_FILE"
echo "BED in: $FEAT_BED_IN"
echo "Base cov file: $BASE_COV_FILE"
echo "BED out: $FEAT_BED_OUT"
echo "nMapUniq: $nMapUniq"
echo "Bedtools options: $S_OPTION"
echo

if [ "$nMapUniq" = "" ]
then
	echo "Did not obtain nMapUniq from stats file. Quitting."
	exit 1
fi


### determine number of columns in input (override with cut by setting in environment)

echo "Getting column number..."
echo

if [[ $FEAT_BED_IN == *".bgz" ]]; then
	N_FEAT_COL=`gunzip -c $FEAT_BED_IN | awk '{print NF; exit}'`
else
	N_FEAT_COL=`awk '{print NF; exit}' $FEAT_BED_IN`
fi


echo "nCols: $N_FEAT_COL"
echo


### run main command

echo "Calculating coverage..."
echo

if [[ $FEAT_BED_IN == *".bgz" ]]; then
	# extract base runs overlapping features and count coverage
	bedtools intersect -nonamecheck -wao $S_OPTION \
	-a <(gunzip -c $FEAT_BED_IN |
	 cut -f 1-$N_FEAT_COL |
	 Rscript $PFX/filter_chroms.R) \
	-b <(gunzip -c $BASE_COV_FILE |
	 awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,".",$5,$4}' |
	 Rscript $PFX/filter_chroms.R) |
	perl $PFX/feature_coverage.pl $nMapUniq $N_FEAT_COL |
	bgzip -c > $FEAT_BED_OUT
	sleep 10
else
	# extract base runs overlapping features and count coverage
	bedtools intersect -nonamecheck -wao $S_OPTION \
	-a <(cat $FEAT_BED_IN |
	 cut -f 1-$N_FEAT_COL |
	 Rscript $PFX/filter_chroms.R) \
	-b <(gunzip -c $BASE_COV_FILE |
	 awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,".",$5,$4}' |
	 Rscript $PFX/filter_chroms.R) |
	perl $PFX/feature_coverage.pl $nMapUniq $N_FEAT_COL > $FEAT_BED_OUT
	sleep 10
fi


echo 
echo "Done!"
