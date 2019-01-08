#!/bin/bash


# Assume both of below files to be n genes x T timepoints
# RAND should be DATAFILE where each gene is randomized by time.


#export DATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/genes.txt
#export NROWS=100     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=insilico_size100_1_1000bootstrap            #protein-codin5
#export NORMALIZATION=0mean
#export LOADREPS=1

#export DATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_3/0mean/reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_3/0mean/reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_3/0mean/genes.txt
#export NROWS=100     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=insilico_size100_3_1000bootstrap            #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1


# WRONG

export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/over-expression-fc/data/KLF15/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps.txt
export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/over-expression-fc/data/KLF15/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps-rand.txt
export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/over-expression-fc/data/KLF15/0mean/genes.txt
export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is only used to determine how to parallelize
export GENES=KLF15-pairwise           #protein-coding
export NORMALIZATION=0mean
export LOADREPS=1

# Just some names for the outputs you use
export DEG=er                   #edge-R
export SAMPLE=norm              #normalized values (or "avg" for averaged values, or "sample" for sampled replicates)


if [ $LOADREPS -eq "1" ]; then
    export REPS=reps;
else
    export REPS=noreps;
fi




# Causal params
export CAUSAL=enet
export LAG=2
export NULL=g


export BOOTSTRAPNUM=1000
#export BOOTSTRAPNUM=100
#export BOOTSTRAPNUM=10

# Output names
export OUTPUTNAME="$GENES-$NORMALIZATION-$REPS-$DEG-$SAMPLE-$CAUSAL-$LAG-$NULL"
export FOLDER=../runs/$OUTPUTNAME                        # Folder made on your directory to scp to cluster.



# Internal processing for the data


export DATANAME=${DATAFILE##*/}  # retain the part after the last slash
export RANDDATANAME=${RANDDATAFILE##*/}
export GENENAME=${GENEFILE##*/}
export HYPERFILE="$CAUSAL/hyperlist_$CAUSAL.p"
export TEST="$(echo $CAUSAL | head -c 1)"


## Job params

# For global null


if [ $NULL == "l" ]; then
    export NPERSCRIPT=1    # Number of rows to run per script.
    export NNODES=1                        # Nodes to run on
    export PPN=4                          # Processors per node.
    export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
    export MEMUPPERBOUND=45000             # In KB
else

    if [ $NULL == "g" ]; then
#        export NPERSCRIPT=15    # Number of rows to run per script.
        export NPERSCRIPT=120    # Number of rows to run per script.
#        export NPERSCRIPT=150    # Number of rows to run per script.
        export NNODES=1                        # Nodes to run on, per script
        export PPN=4                          # Processors per node, per script
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





# Just some echo'ing so you can confirm we have the right values

echo FOLDER is $FOLDER
echo DATAFILE is $DATAFILE
echo RANDDATAFILE is $RANDDATAFILE
echo DATANAME is $DATANAME
echo RANDDATANAME is $RANDDATANAME
echo TEST is $TEST
echo OUTPUTNAME is $OUTPUTNAME

#
echo NROWS is $NROWS
echo SCRIPTNUM is $SCRIPTNUM
echo PARALLELNUM is $PARALLELNUM
echo TIMEUPPERBOUND is $TIMEUPPERBOUND
echo MEMUPPERBOUND is $MEMUPPERBOUND
