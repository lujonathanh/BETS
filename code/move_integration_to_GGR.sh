#!/usr/bin/env bash

source ./package_params_cpipeline.sh

#topdir=/Users/jlu96/v-causal-snps/code/della/test-over-expression-dir
topdir=/tigress/BEE/RNAseq/Data/Networks/GGR/integration

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

ls -d run_*-fdr > run_folders.txt

while read runf; do
    pref="run_"
    folder=$topdir/$NORM"_"$NULL-null"_"${runf#${pref}}

    mkdir $folder
    mkdir $folder/run

    mkdir $folder/run/$GENES-$CAUSAL-$LAG
    mkdir $folder/networks
    echo "Moving " $runf "to" $folder/run/$GENES-$CAUSAL-$LAG
    cp -r $runf/. $folder/run/$GENES-$CAUSAL-$LAG
    echo "Moving " $runf " union networks to " $folder/networks
    cp -r $runf/networks/*union* $folder/networks
done < run_folders.txt