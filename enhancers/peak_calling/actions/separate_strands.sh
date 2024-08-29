#!/bin/bash

### description

# script written by Brian Magnuson; script modified by Ariel McShane (20240819)

SAMPLE=$1
STRAND=$2

GENOME_BAM=$SAMPLE_DIR/$SAMPLE.$BAM_SUFFIX.bam
STRAND_BAM=$PROJECT_DIR/$SAMPLE.$BAM_SUFFIX.$STRAND.bam

function separate_strands {

	if [[ "$IS_PAIRED" == "TRUE" ]]
	then
		echo "$SAMPLES has paired-end reads."
		if [[ "$STRAND" == "plus" ]]
		then
			F="(not secondary_alignment and proper_pair) and ((first_of_pair and reverse_strand and template_length > -1000 and template_length < -150) or (second_of_pair and not reverse_strand and template_length < 1000 and template_length > 150))"
		else
			# strand is minus
			F="(not secondary_alignment and proper_pair) and ((first_of_pair and not reverse_strand and template_length < 1000 and template_length > 150) or (second_of_pair and reverse_strand and template_length > -1000 and template_length < -150)) "
		fi
	else
		echo "$SAMPLES has single-end reads."
		if [[ "$STRAND" == "plus" ]]
		then
			F="(not (secondary_alignment or unmapped) and reverse_strand)"
		else
			# strand is minus
			F="(not (secondary_alignment or unmapped) and not reverse_strand)"
		fi
	fi
	sambamba view \
		-F "$F" \
		-f bam \
		-h \
		-t 16 \
		$GENOME_BAM \
	> $STRAND_BAM
	
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

#######################################
# sep strand into new bam files
separate_strands
#######################################

# ### unload conda env
# conda deactivate
# # reactivate old env (if applicable)
# if [[ "$OLD_ENV" != "" ]]
# then
	# conda activate $OLD_ENV
# fi
