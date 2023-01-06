#!/bin/bash
set -e

source ./package_params_cpipeline.sh

python3 prep_jobs_rand_cv.py -d $DATANAME -d2 $RANDDATANAME -lr $LOADREPS -o $OUTPUTNAME  -hlf $HYPERFILE -t $TEST -tn $CAUSAL -sn $SCRIPTNUM -p $PARALLELNUM -l $LAG -n $NULL

echo
echo Done prepping jobs