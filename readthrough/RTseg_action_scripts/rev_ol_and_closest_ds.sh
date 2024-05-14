#!/bin/bash

# This script uses bedtools to find RT-segment gene overlaps on the opposite strand and to identify the
# distance of a RT segment to the nearest downstream gene on the same strand for classification.

#q    require $SAMPLE $GB_RAM $GENE_ANNOT
#q    require $OUTDIR

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb


# unzip count_bed output
gunzip -c $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed.bgz > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed

# get reverse overlaps
echo "count reverse strand gene overlaps"
bedtools intersect -S -wao -a $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed -b $GENE_ANNOT > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.reverseOL.bed

# get closest ds genes
echo "find distance to closest downstream gene"
bedtools closest -s -D a -iu -a $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.fixedRPKM.bed -b $GENE_ANNOT > $OUTDIR/$SAMPLE.RTsegments.TESintersect.segmentRPKM.closestDS.bed
