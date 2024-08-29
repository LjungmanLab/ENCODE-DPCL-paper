#!/bin/bash

### description

# script written by Brian Magnuson using modified commands from Ramil Nurtdinov; script modified by Ariel McShane (20240819)


STRAND="$1"; shift
SAMPLES=("$@")

unset STRAND_BAMS
STRAND_BAMS=$( \
	for SAMPLE in ${SAMPLES[@]}
	do
	
		echo -n " $PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.$STRAND.bam"
		
	done)


if [[ "$CALL_SUMMITS" == "TRUE" ]]
then
	PEAKS_OUTDIR=$PROJECT_DIR/$UV_GROUP.$STRAND.summits
else
	PEAKS_OUTDIR=$PROJECT_DIR/$UV_GROUP.$STRAND
fi


if [[ "$GENOME" == "hg38" ]] || [[ "$GENOME" == "hg19" ]] || [[ "$GENOME" == "hg18" ]]
then
	GENOME_MACS2=hs
elif [[ "$GENOME" == "mm9" ]] || [[ "$GENOME" == "mm10" ]]
then
	GENOME_MACS2=mm
fi

function call_peaks_replicates {
	
	if [[ "$CALL_SUMMITS" == "TRUE" ]]
	then
		SUMMITS_CMD="--call-summits"
	else
		SUMMITS_CMD=""	
	fi

	echo "calling peaks for $UV_GROUP $STRAND..."
	macs2 callpeak \
		-t $STRAND_BAMS \
		-f BAMPE \
		-g $GENOME_MACS2 \
		--outdir $PEAKS_OUTDIR \
		-n $UV_GROUP.$STRAND \
		-p 0.01 \
		--SPMR \
		--bdg \
		--keep-dup 1 \
		$SUMMITS_CMD
}


echo
echo "---- checks ----"
echo $STRAND_BAMS
echo $CALL_SUMMITS
echo $GENOME_MACS2
echo $PEAKS_OUTDIR
echo $UV_GROUP.$STRAND
echo "----------------"
echo
	
	
# ### load conda env

# OLD_ENV=""

# echo "conda: $CONDA_DEFAULT_ENV"

# if [ -f "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh" ]
# then
	# . "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh"
# fi


# if [[ "$CONDA_DEFAULT_ENV" != "mac2" ]]
# then
	# if [[ "$CONDA_DEFAULT_ENV" == "base" ]]
	# then
		# conda deactivate
	# elif [[ "$CONDA_DEFAULT_ENV" != "" ]]
	# then
		# OLD_ENV=$CONDA_DEFAULT_ENV
		# conda deactivate
	# fi
	# conda activate macs2
# fi

# print software versions
sambamba --version
macs2 --version

#######################################
# call peaks using replicates
call_peaks_replicates
#######################################

# ### unload conda env
# conda deactivate
# # reactivate old env (if applicable)
# if [[ "$OLD_ENV" != "" ]]
# then
	# conda activate $OLD_ENV
# fi
