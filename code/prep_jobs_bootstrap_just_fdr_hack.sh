#!/bin/bash

source ./package_params_cpipeline.sh

python2 prep_jobs_bootstrap_just_fdr_hack.py -d $DATANAME -d2 $RANDDATANAME -lr $LOADREPS -o $OUTPUTNAME -t $TEST -tn $CAUSAL -sn $SCRIPTNUM -p $PARALLELNUM -l $LAG -n $NULL -bn $BOOTSTRAPNUM

echo
echo Done prepping jobs