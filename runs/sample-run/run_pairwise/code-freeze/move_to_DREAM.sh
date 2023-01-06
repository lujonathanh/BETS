#!/usr/bin/env bash

source ./package_params_cpipeline.sh



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

topdir=/tigress/jhlu/DREAM/results


export runf="run_l-fdr"

folder=$topdir/$GENES

mkdir $folder
mkdir $folder/run
mkdir $folder/run/$NORM-$CAUSAL-$LAG
mkdir $folder/networks
echo "Moving " $runf "to" $folder/run/$NORM-$CAUSAL-$LAG
cp -r $runf/. $folder/run/$NORM-$CAUSAL-$LAG
echo "Moving " $runf " union networks to " $folder/networks
cp -r $runf/networks/*union* $folder/networks

mkdir $folder/bootstrap-results
mkdir $folder/bootstrap-results-fdr-0.05-effect

echo "Moving " $runf " bootstrap integrated networks to " $folder/bootstrap-results
cp -r $runf/bootstrap-results/*union* $folder/bootstrap-results
echo "Moving " $runf " fdr-thresholded bootstrap integrated networks to " $folder/bootstrap-results-fdr-0.05-effect
cp -r $runf/bootstrap-results-fdr-0.05-effect/*union* $folder/bootstrap-results-fdr-0.05-effect