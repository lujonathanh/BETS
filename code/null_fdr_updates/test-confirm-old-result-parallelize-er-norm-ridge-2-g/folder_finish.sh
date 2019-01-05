#!/usr/bin/env bash

source ./package_params_cpipeline.sh

mkdir finish
mkdir finish/hyper
mkdir finish/fit

cp *summary* finish
cp -r fdr-0.05* finish
cp -r plots finish
cp -r hyper/*_* finish/hyper
cp -r fit/$OUTPUTNAME"_fit"* finish/fit

echo "Moving main files to folder: finish/"