#!/bin/bash

source ./package_params_cpipeline.sh

python3 prep_jobs_bootstrap.py -d $DATANAME -d2 $RANDDATANAME -lr $LOADREPS -o $OUTPUTNAME -t $TEST -tn $CAUSAL -sn $SCRIPTNUM -p $PARALLELNUM -l $LAG -n $NULL -bn $BOOTSTRAPNUM

echo
echo Done prepping jobs