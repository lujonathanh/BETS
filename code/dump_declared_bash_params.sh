#!/bin/bash


export SCRIPT="package_params_cpipeline.sh"
export OUTPUT="params.csv"
echo "Dumping from " $SCRIPT " to " $OUTPUT

source ./package_params_cpipeline.sh

python2 get_declared_bash_params.py -f $SCRIPT -o tmp

printenv > tmp-envs.txt

python2 dump_bash_params.py -e tmp-envs.txt -p tmp -o $OUTPUT

rm tmp
rm tmp-envs.txt