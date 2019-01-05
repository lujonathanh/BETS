#!/bin/bash

echo Make sure to update package_params_cpipeline.sh first!

source ./package_params_cpipeline.sh

mkdir $FOLDER

cp $DATAFILE $FOLDER
cp $RANDDATAFILE $FOLDER
cp -r $CAUSAL $FOLDER

while read file; do
    echo Copy $file to $FOLDER
    cp $file $FOLDER
done <package_required_files_cpipeline.txt


if [ "$LOADREPS" != "0" ];
    then
        export DIR=$(dirname "${DATAFILE}")

        while read file; do
            echo Copy $DIR/$file to $FOLDER
            cp $DIR/$file $FOLDER
        done <$DATAFILE

        export RANDDIR=$(dirname "${RANDDATAFILE}")

        while read file; do
            echo Copy $RANDDIR/$file to $FOLDER
            cp $RANDDIR/$file $FOLDER
        done <$RANDDATAFILE
fi


