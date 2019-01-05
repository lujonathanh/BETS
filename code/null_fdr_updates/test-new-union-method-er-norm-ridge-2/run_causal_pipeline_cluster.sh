#!/usr/bin/env bash

source ./package_params_cpipeline.sh
export script="time python causal_pipeline.py -df $DATANAME -rdf $RANDDATANAME -o $OUTPUTNAME -hl $HYPERFILE -t $TEST -l $LAG -sb $STRATIFY"

echo Submitting Parallel Script $script
clusterize -l $TIMEUPPERBOUND:00 -m $MEMUPPERBOUND -n $NNODES -p $PPN  -c "time $script & wait"
