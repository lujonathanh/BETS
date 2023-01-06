#!/usr/bin/env bash

source ./package_params_cpipeline.sh

#topdir=/Users/jlu96/v-causal-snps/code/della/test-over-expression-dir
topdir="/tigress/BEE/RNAseq/old/Data/Networks/GGR/pairwise"

mkdir $topdir

if [ "$NORMALIZATION" = "0mean1var" ];
then
    export NORM="0mean-1var"
else
    if [ "$NORMALIZATION" = "0mean" ];
    then
        export NORM="0mean-unnormalized"
    else
        echo "Error: NORMALIZATION $NORMALIZATION is not 0mean or 0mean1var"
        exit 1
    fi
fi


folder=$topdir/$NORM"_pairwise"

echo "Moving run_pairwise to" $folder/run/$GENES-pairwise-$LAG
cp -r run_pairwise/. $folder/run/$GENES-pairwise-$LAG
