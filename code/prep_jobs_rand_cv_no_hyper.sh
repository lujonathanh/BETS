#!/bin/bash

source ./package_params_cpipeline.sh

python2 prep_jobs_rand_cv.py -d $DATANAME -d2 $RANDDATANAME -lr $LOADREPS -o $OUTPUTNAME  -hlf $HYPERFILE -t $TEST -tn $CAUSAL -sn $SCRIPTNUM -p $PARALLELNUM -l $LAG -n $NULL -cv 0

echo
echo Done prepping jobs