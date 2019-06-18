#!/bin/bash

source ./package_params_cpipeline.sh

python2 prep_jobs_pairwise.py -d $DATANAME -lr $LOADREPS -o $OUTPUTNAME -sn $SCRIPTNUM -p $PARALLELNUM -l $LAG

echo
echo Done prepping jobs