#!/bin/bash

### call peaks for BruUV-seq samples and up to 2 replicates


### assign batch headers as needed

#SBATCH --job-name=bruuv_peak_calling		# can either leave this job as generic or make specific (specific names will help find /log file more easily)
#SBATCH --mail-user=mcshanea@umich.edu		# add user email address or comment out if update communications not desired
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --nodes=1							# will use this number of nodes per job submitted, thus 1 is likely sufficient for most applications
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --account=ljungman0
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.log


### General variables (***will simplify this for sharing purposes since all filenames will have same components)
# ***simplifying naming conventions/variables will need to be propigated through all scripts

CONDA_ENV=/nfs/turbo/radonc-ljungman-turbo/envs/macs2							# path to conda environment

# set directories
WORKING_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/peak_calling/working
ACTIONS_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/peak_calling/actions
SAMPLE_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/peak_calling/sample		# location of .bam files
PROJECT_DIR=$WORKING_DIR/peak_calling_test

# stranded bam file names will be formated as $SAMPLE.$BAM_SUFFIX.$STRAND.bam
BAM_SUFFIX="hg38.gencode_29.spikeins_1.genome"

# other variables
GENOME=hg38						# needed to determine -g option in callpeak command
CHROM_SIZES=/nfs/turbo/radonc-ljungman-turbo/mcshanea/ENCODE/genome_cov/chrom_sizes/hg38.spikeins_1.EBV_ERCC_phiX.UCSC.genome.chromSizes.txt
IS_PAIRED=TRUE					# TRUE if bam files are paired-end


### create project directory

mkdir -p $PROJECT_DIR


# ***each set of samples was run like this in parallel... figure out best way to set this up in shell
SAMPLES1=("HCT1160h4001a" "HCT1160h4002a")
UV_GROUP1=HCT116UV

. $ACTIONS_DIR/bruuvseq_peak_calling_main_workflow.sh $UV_GROUP1 "${SAMPLES1[@]}"



