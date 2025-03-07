#!/bin/bash
DATE="18APR2023"

CURR_DIR=`pwd`
BT_STAT=("wo_breakthrough" "w_breakthrough")

INF_STAT=("non_inf" "inf")

B_ADJ="D85D147D237"

for bt in ${BT_STAT[*]}; do

for inf in ${INF_STAT[*]}; do

echo $inf
cd $CURR_DIR
cd ./figures/${DATE}/stage1/${inf}/landscapes/${B_ADJ}_b+o_adj/${bt}/optimization_1

PNG_LIST=`ls -lart *.png | awk '{print $9}'`

echo "Cropping landscapes"

for f in $PNG_LIST; do
 #   convert ${f} -crop 2200x1200+500+600 ${f}
    convert ${f} -gravity center -crop 1950x1100-140+150 ${f}
done

done

done
