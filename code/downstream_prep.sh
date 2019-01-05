#!/usr/bin/env bash

# The goal here is to aggregate the downstream files for networks with 0.05 for ease of porting. Includes the hyperparameters and fitted networks

source ./package_params_cpipeline.sh

FDRS=("g" "l")
FDRALIASES=("none" "effect")
FDRTHRESHOLD=0.05

for i in `seq 0 1`;
do
    FDR=${FDRS[i]}
    FDRALIAS=${FDRALIASES[i]}
    echo $FDR
    echo $FDRALIAS

    OUTDIR=run_$FDR-fdr


    if [ -z "$(ls -A fdr-$FDRTHRESHOLD-$FDRALIAS)" ]; then
        echo fdr-$FDRTHRESHOLD-$FDRALIAS is empty
    else

        mkdir $OUTDIR

        echo "Moving main files to folder: $OUTDIR"

        echo "### DATA ###"

        mkdir $OUTDIR/data
        cp $DATANAME $OUTDIR/data
        cp $RANDDATANAME $OUTDIR/data
        cp $GENENAME $OUTDIR/data

        touch $OUTDIR/data/files.csv
        echo DATANAME,$DATANAME >> $OUTDIR/data/files.csv
        echo RANDDATANAME,$RANDDATANAME >> $OUTDIR/data/files.csv

        if [ "$LOADREPS" != "0" ];
        then
            while read file; do
                cp $file $OUTDIR/data
            done <$DATANAME

            while read file; do
                cp $file $OUTDIR/data
            done <$RANDDATANAME
        fi

        echo "### HYPER ###"
        mkdir $OUTDIR/hyper
        cp -r hyper/*_* $OUTDIR/hyper

        echo "### FIT ###"
        mkdir $OUTDIR/fit
        cp -r fit/$OUTPUTNAME"_fit"* $OUTDIR/fit
        cp fit_all_summary_normal.txt $OUTDIR/fit
        cp fit_all_summary_fdr-$FDRALIAS.txt $OUTDIR/fit/fit_all_summary_$FDR-fdr.txt

        echo "### NETWORKS ###"
        mkdir $OUTDIR/networks
        ls fdr-$FDRTHRESHOLD-$FDRALIAS > tmp
        while read file; do
            cp fdr-$FDRTHRESHOLD-$FDRALIAS/$file $OUTDIR/networks
            #cp -r fdr-$FDRTHRESHOLD-$FDRALIAS/$file $OUTDIR/networks/$(python della_convert_filename.py -f $file)
        done < tmp

        rm tmp

        echo "### RAW NETWORK ### "
        mkdir $OUTDIR/raw-network
        cp *-union-network.txt $OUTDIR/raw-network

        echo "### PLOTS ###"
        cp -r plots $OUTDIR

        echo "### TIMING ###"
        cp -r timing $OUTDIR

        echo "### CODE FREEZE ###"
        mkdir $OUTDIR/code-freeze
        while read file; do
            cp $file $OUTDIR/code-freeze
        done < package_required_files_cpipeline.txt


        echo "### PARAMS ###"
        cp package_params_cpipeline.sh $OUTDIR
        cp params.csv $OUTDIR




        echo "### BOOTSTRAP ###"

        BOOTSTRAPFOLDER=bootstrap/bootstrap-results

        if [ -z "$(ls -A $BOOTSTRAPFOLDER)" ]; then
           echo $BOOTSTRAPFOLDER is empty
        else
           cp -r $BOOTSTRAPFOLDER $OUTDIR
           echo $BOOTSTRAPFOLDER copied to $OUTDIR
        fi

        BOOTSTRAPFDRFOLDER=bootstrap/bootstrap-results-fdr-$FDRTHRESHOLD-$FDRALIAS

        if [ -z "$(ls -A $BOOTSTRAPFDRFOLDER)" ]; then
            echo $BOOTSTRAPFDRFOLDER is empty
        else
          cp -r $BOOTSTRAPFDRFOLDER $OUTDIR
          echo $BOOTSTRAPFDRFOLDER copied to $OUTDIR
        fi
    fi

done




# Where you copy the pairwise part.

OUTDIR=run_pairwise

mkdir $OUTDIR

echo "Moving main files to folder: $OUTDIR"

echo "### DATA ###"

mkdir $OUTDIR/data
cp $DATANAME $OUTDIR/data
cp $GENENAME $OUTDIR/data

touch $OUTDIR/data/files.csv
echo DATANAME,$DATANAME >> $OUTDIR/data/files.csv

if [ "$LOADREPS" != "0" ];
then
    while read file; do
        cp $file $OUTDIR/data
    done <$DATANAME
fi

echo "### PAIRWISE ###"
PAIRWISEFOLDER=pairwise/pairwise-results

if [ -z "$(ls -A $PAIRWISEFOLDER)" ]; then
   echo $PAIRWISEFOLDER is empty
else
   cp -r $PAIRWISEFOLDER $OUTDIR
   echo $PAIRWISEFOLDER copied to $OUTDIR
fi

echo "### TIMING ###"
cp -r timing $OUTDIR

echo "### CODE FREEZE ###"
mkdir $OUTDIR/code-freeze
while read file; do
    cp $file $OUTDIR/code-freeze
done < package_required_files_cpipeline.txt


echo "### PARAMS ###"
cp package_params_cpipeline.sh $OUTDIR
cp params.csv $OUTDIR