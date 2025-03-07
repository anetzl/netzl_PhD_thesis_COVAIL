#!/bin/bash

DATE="18APR2023"

CURR_DIR=`pwd`

BT_STAT=("w_breakthrough")

#INF_STAT=("non_inf" "inf")

inf=${1}

PRE=${2}

# B_ADJ="D85D147D237"

B_ADJ="wo_2dose"

for bt in ${BT_STAT[*]}; do

#for inf in ${INF_STAT[*]}; do

echo $inf
cd $CURR_DIR
# cd ./figures/${DATE}/stage1/${inf}/landscapes/${B_ADJ}_b+o_adj/${bt}/optimization_1


cd ./figures/${DATE}/stage1/${inf}/landscapes_imprinting/${B_ADJ}_b+o_adj/${bt}/optimization_1/${PRE}

HTML_LIST=`ls -lart *.html | awk '{print $9}'`

echo "Capturing screenshots of html landscapes"

for f in $HTML_LIST; do
    new_name=${f%.*}.png
    echo $new_name
    open $f
    sleep 1
   # screencapture -x -w $new_name
    screencapture -x -m $new_name
#done

done

done
