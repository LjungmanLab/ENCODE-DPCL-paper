#!/bin/bash

### description

# script written by Brian Magnuson using modified commands from Ramil Nurtdinov; script modified by Ariel McShane (20240819)

SAMPLE=$1
STRAND=$2

STRAND_BAM=$PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.$STRAND.bam


if [[ "$CALL_SUMMITS" == "TRUE" ]]
then
	PEAKS_OUTDIR=$PROJECT_DIR/$SAMPLE.$STRAND.summits
else
	PEAKS_OUTDIR=$PROJECT_DIR/$SAMPLE.$STRAND
fi


if [[ "$GENOME" == "hg38" ]] || [[ "$GENOME" == "hg19" ]] || [[ "$GENOME" == "hg18" ]]
then
	GENOME_MACS2=hs
elif [[ "$GENOME" == "mm9" ]] || [[ "$GENOME" == "mm10" ]]
then
	GENOME_MACS2=mm
fi

function call_peaks {
	
	SUMMITS="$1"
	
	if [[ "$SUMMITS" == "TRUE" ]]
	then
		SUMMITS_CMD="--call-summits"
	else
		SUMMITS_CMD=""	
	fi

	echo "calling peaks for $SAMPLE $STRAND..."
	macs2 callpeak \
		-t $STRAND_BAM \
		-f BAMPE \
		-g $GENOME_MACS2 \
		--outdir $PEAKS_OUTDIR \
		-n $SAMPLE.$STRAND \
		-p 0.01 \
		--SPMR \
		--bdg \
		--keep-dup 1 \
		$SUMMITS_CMD
}
		
# ### load conda env
# echo "conda: $CONDA_DEFAULT_ENV"

# if [ -f "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh" ]
# then
	# . "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh"
# fi


# OLD_ENV=""
# if [[ "$CONDA_DEFAULT_ENV" != "macs2" ]]
# then
	# if [[ "$CONDA_DEFAULT_ENV" == "base" ]]
	# then
		# conda deactivate
	# elif [[ "$CONDA_DEFAULT_ENV" != "" ]]
	# then
		# OLD_ENV=$CONDA_DEFAULT_ENV
		# conda deactivate
	# # else
		# # if [ -f "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh" ]
		# # then
			# # . "/garage/wilsonte_lab/bin/anaconda3/etc/profile.d/conda.sh"
		# # fi
	# fi
	# conda activate macs2
# fi

# print software versions
sambamba --version
macs2 --version

#######################################
# call peaks on individual sample
call_peaks "$CALL_SUMMITS" 
#######################################


# ### unload conda env
# conda deactivate
# # reactivate old env (if applicable)
# if [[ "$OLD_ENV" != "" ]]
# then
	# conda activate $OLD_ENV
# fi
