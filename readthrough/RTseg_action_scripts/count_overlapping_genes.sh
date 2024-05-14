#!/bin/bash

# This script uses bedtools to create the overlapping gene file for rt_segment evaluation. The resulting
# file is then zipped for use in count_bed.

#q    require $SAMPLE $GB_RAM $GENE_ANNOT
#q    require $IN_FILE $OUT_FILE

#$    -N  $STEP_NAME\_$SAMPLE
#$    -wd $OUTDIR
#$    -l  vf=$GB_RAM\G

#PBS  -N  $STEP_NAME\_$SAMPLE
#PBS  -d  $OUTDIR
#PBS  -l  mem=$GB_RAM\gb

echo "count all gene overlaps"
bedtools intersect -s -wao -a $IN_FILE -b $GENE_ANNOT > $OUT_FILE

awk 'BEGIN{ FS = OFS = "\t" } { print $0, "X" }' $OUT_FILE > tmp && mv tmp $OUT_FILE

gzip -c $OUT_FILE > $OUT_FILE\.bgz
