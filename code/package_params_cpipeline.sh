#!/bin/bash





# For every run, you must set the values under "DATA INPUT"
# and "NAMES".



####################################
##           DATA INPUT           ##
####################################


export LOADREPS=1
# Set to 1 if there are replicates. Each replicate should be a complete set of gene expression timeseries over the same number of timepoints.
# Set to 0 if there is only one data set.


export DATAFILE=../data/DREAM/insilico_size100_1/0mean/reps.txt
# If there are replicates, this should go to
# If no replicates, this should go to a file that has n genes and T timepoints. It should have the first column be the genenames, and the columns after be the ordered timepoints.


export RANDDATAFILE=../data/DREAM/insilico_size100_1/0mean/reps-rand.txt
# Same as DATAFILE, but the files should be the original files where each gene has been randomized over time.


export GENEFILE=../data/DREAM/insilico_size100_1/0mean/genes.txt
# List of gene names corresponding to each row of the data files.


export NROWS=100
# Number of rows (i.e. # genes) in your datafile.
# This is only used to determine the number of scripts to make when parallelizing.




####################################
##             NAMES              ##
####################################



export GENES=insilico_size100_1-test-namedellamodify
# Give a name of the dataset.


export NORMALIZATION=0mean
# Give a name to the normalization, e.g. 0mean for zero-mean centered data.
# BETS does NOT do normalization for you! We recommend zero-mean centered data.


export DEG=er
# Just a suffix that is used (TODO: REMOVE)


export SAMPLE=norm
# Just a name for the SAMPLE (TODO: REMOVE)




####################################
##           PARAMETERS           ##
####################################
# Parameters of the causal network inference method.


export CAUSAL=lasso
# The type of vector autoregression model.
# Must be one of {ridge, lasso, enet}
# BETS sets CAUSAL=enet


export LAG=2
# The lag of vector autoregression model.
# Must be a positive integer that is less than the number of timepoints you have.
# BETS sets LAG=2


export BOOTSTRAPNUM=1000
# Number of bootstrap samples to use for stability selection
# Set to a positive integer, e.g. 100 or 1000
# BETS sets BOOTSTRAPNUM=1000
# If your time is limited, BOOTSTRAPNUM=100 should be fine



############ NO NEED TO MODIFY ANY CODE FROM BELOW #################
# Internal processing for the data


export OUTPUTNAME="$GENES-$NORMALIZATION-$REPS-$DEG-$SAMPLE-$CAUSAL-$LAG-$NULL"
# Name of the folder


export FOLDER=../runs/$OUTPUTNAME
# Where the folder will be copied to. We recommend keeping all the runs in a single place.


export NULL=g  # DO NOT CHANGE
# The null thresholding method, where g is broad null and l is the narrow null
# Must be one of {l, g}.
# BETS uses g.


###################### JOB PARAMETERS ###############################

if [ $NULL == "g" ]; then
export NPERSCRIPT=120                  # Number of rows to run per script.
export NNODES=1                        # Nodes to run on, per script
export PPN=4                          # Processors per node, per script
export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
export MEMUPPERBOUND=45000             # In KB

else

if [ $NULL == "l" ]; then
export NPERSCRIPT=1     # Number of rows to run per script.
export NNODES=1                        # Nodes to run on
export PPN=4                          # Processors per node.
export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
export MEMUPPERBOUND=45000             # In KB
else

echo "NULL is not 'g' or 'l'. Error!!!"
exit 1

fi
fi

if [ "$NPERSCRIPT" -eq "1" ];
    then export SCRIPTNUM=$NROWS;
    else export SCRIPTNUM=$(expr $NROWS / $NPERSCRIPT + 1);
fi

export PARALLELNUM=$(($NNODES * $PPN)) # (Don't change) # of scripts to run per job, one per processor
###############################################################


# Folder to copy

export DATANAME=${DATAFILE##*/}  # retain the part after the last slash
export RANDDATANAME=${RANDDATAFILE##*/}
export GENENAME=${GENEFILE##*/}
export HYPERFILE="$CAUSAL/hyperlist_$CAUSAL.p"
export TEST="$(echo $CAUSAL | head -c 1)"


if [ $LOADREPS -eq "1" ]; then
export REPS=reps;
else
export REPS=noreps;
fi


# Just some echo'ing so you can confirm we have the right values

echo FOLDER is $FOLDER
echo DATAFILE is $DATAFILE
echo RANDDATAFILE is $RANDDATAFILE
echo CAUSAL is $CAUSAL
echo OUTPUTNAME is $OUTPUTNAME

#
echo NROWS is $NROWS
echo SCRIPTNUM is $SCRIPTNUM
echo PARALLELNUM is $PARALLELNUM
echo TIMEUPPERBOUND is $TIMEUPPERBOUND
echo MEMUPPERBOUND is $MEMUPPERBOUND