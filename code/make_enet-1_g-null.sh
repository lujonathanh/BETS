#!/usr/bin/env bash

cp package_params_cpipeline.sh tmp

sed 's/export CAUSAL=.*/export CAUSAL=enet/g; s/export LAG=.*/export LAG=1/g; s/export NULL=.*/export NULL=g/g;' tmp > package_params_cpipeline.sh

source ./package_params_cpipeline.sh
./package_for_cluster_cpipeline.sh


cp tmp package_params_cpipeline.sh
rm tmp