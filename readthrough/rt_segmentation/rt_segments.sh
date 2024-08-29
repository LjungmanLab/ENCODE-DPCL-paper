#!/bin/bash

### Perform readthrough segments analysis as executed in McShane et al. Primary input files are
### Bru-seq 0h 250bp segments files (.bed). These and accessory files required for running this
### workflow are hosted on AWS. These files must be downloaded to the working directory.


### assign batch headers as needed

#SBATCH --job-name=DPCL_readthrough			# can either leave this job as generic or make specific (specific names will help find /log file more easily)
#SBATCH --mail-user=mcshanea@umich.edu		# add user email address or comment out if update communications not desired
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --nodes=1							# will use this number of nodes per job submitted, thus 1 is likely sufficient for most applications
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
#SBATCH --account=ljungman0
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.log
#SBATCH --array=0-1							# be sure to set array to the number of samples (0-N)


### General variables (***will simplify this for sharing purposes since all filenames will have same components)
### ***Add option for keeping or removing intermediate files

CONDA_ENV=/nfs/turbo/radonc-ljungman-turbo/envs/p39							# path to conda environment

PROJECT=readthru_test
WORKING_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/readthru/working		# put all segments, annotation, and counts/rpkm files here
ACTIONS_DIR=/nfs/turbo/radonc-ljungman-turbo/mcshanea/code_test/readthru/actions

OUT_DIR=$WORKING_DIR/$PROJECT

BRU_TP='0h'
BIN_SIZE='250'

SEG_SUFFIX=hg38.gencode_29.spikeins_1.bin_$BIN_SIZE.segments.actual_RPKM.bed.bgz
GENE_ANNOT=$WORKING_DIR/gencode.v29.basic.annotation_genes_final_sorted.bed
TES_ANNOT=$WORKING_DIR/gencode.v29.basic.annotation_tes_final_sorted.bed
COUNTS_0h=$WORKING_DIR/ALLcounts_ENCODE16CL_0h_genes.bed
RPKM_0h=$WORKING_DIR/ALLrpkm_ENCODE16CL_0h_genes.bed

# For general use with no BruChase-seq 6h samples, set RPKM_6h=NONE script will use 0h rpkm file for classification
RPKM_6h=$WORKING_DIR/ALLrpkm_ENCODE16CL_6h_exons.bed
# RPKM_6h=NONE

# Remove intermediate files?
REMOVE_FILES=TRUE


### Samples array

# SAMPLES=("A6730h4001a" "A6730h4002a" "Caco20h4001a" "Caco20h4002a"
		 # "Calu30h4001a" "Calu30h4002a" "GM128780h4001a" "GM128780h4002a"
		 # "HCT1160h4001a" "HCT1160h4002a" "HepG20h4005a" "HepG20h4006a"
		 # "HMEC0h4001a" "HMEC0h4002a" "HUVEC0h4001a" "HUVEC0h4002a"
		 # "IMR900h4003a" "IMR900h4004a" "K5620h4004a" "K5620h4005a"
		 # "MCF10A0h4001a" "MCF10A0h4002a" "MCF70h4001a" "MCF70h4002a"
		 # "OCILY70h4001a" "OCILY70h4002a" "panc10h4001a" "panc10h4002a"
		 # "PC30h4001a" "PC30h4002a" "PC90h4001a" "PC90h4002a"
		 # )

SAMPLES=("HCT1160h4001a"
		 "HCT1160h4002a")
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}		# gets sample name from array based on slurm array index


### Run main workflow

. $ACTIONS_DIR/rt_main_workflow.sh $SAMPLE

