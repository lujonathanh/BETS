#!/bin/bash

source ./package_params_cpipeline.sh
time python causal_pipeline.py -df $DATANAME -rdf $RANDDATANAME -o $OUTPUTNAME -hl $HYPERFILE -t $TEST -l $LAG -sb $STRATIFY