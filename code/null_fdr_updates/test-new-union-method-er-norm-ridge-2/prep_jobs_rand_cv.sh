#!/bin/bash

source ./package_params_cpipeline.sh

python prep_jobs_rand_cv.py -d $DATANAME -d2 $RANDDATANAME -lr $LOADREPS -o $OUTPUTNAME  -hlf $HYPERFILE -t $TEST -tn $CAUSAL -n $SCRIPTNUM -p $PARALLELNUM -l $LAG -sb $STRATIFY

echo
echo Done prepping jobs