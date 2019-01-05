#!/bin/bash


# Assume both of below files to be n genes x T timepoints
# RAND should be DATAFILE where each gene is randomized by time.

#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean1var/small_er-rep3-2cutoff-0mean-1var-100.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean1var/small_er-rep3-2cutoff-0mean-1var-100-rand.txt
#export NROWS=100
#export GENES=test-new-fdr
#export LOADREPS=0


#
export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean1var/small_er-reg-cutoff-2-plus-GR-0mean1var-reps.txt
export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean1var/small_er-reg-cutoff-2-plus-GR-0mean1var-reps-rand.txt
export NROWS=100
export GENES=test-confirm-old-result
export LOADREPS=1


#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean1var/edgeR-reg-cutoff-2-plus-GR-0mean1var-reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean1var/edgeR-reg-cutoff-2-plus-GR-0mean1var-reps-rand.txt
#export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=prot2TPM-0mean1var-reps       #protein-coding
#export LOADREPS=1



#
#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps-rand.txt
#export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=prot2TPM-0mean-reps       #protein-coding
#export LOADREPS=1


# Just some names for the outputs you use
export DEG=er                   #edge-R
export SAMPLE=norm              #normalized values (or "avg" for averaged values, or "sample" for sampled replicates)


# Causal params
export CAUSAL=ridge       # Name of causal folder
export HYPERFILE="$CAUSAL/hyperlist_$CAUSAL.p"
export TEST=r                  # THIS SHOULD CORRESPOND TO CAUSAL
export LAG=2

# Null, FDR params
export NULL=g               # l for local, g for global
# no longer needed: export STRATIFY=n      # Set to "e" to stratify by effect when doing FDR. Else to "n" to not stratify at all. See README for description.


# Output names
export OUTPUTNAME="$GENES-$DEG-$SAMPLE-$CAUSAL-$LAG-$NULL"
export FOLDER=della/$OUTPUTNAME                        # Folder made on your directory to scp to cluster.




export DATANAME=${DATAFILE##*/}  # retain the part after the last slash
export RANDDATANAME=${RANDDATAFILE##*/}


## Job params
#export NPERSCRIPT=20    # Number of rows to run per script.
#export NNODES=4                        # Nodes to run on
#export PPN=4                          # Processors per node.
#export TIMEUPPERBOUND=30              # In minutes. For pairwise, 60 is usually good, for enet, up to 180 needed.
#export MEMUPPERBOUND=45000             # In KB

export NPERSCRIPT=150    # Number of rows to run per script.
export NNODES=1                        # Nodes to run on
export PPN=4                          # Processors per node.
export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
export MEMUPPERBOUND=45000             # In KB


# For local null: to ensure fast enough
#export NPERSCRIPT=50    # Number of rows to run per script.
#export NNODES=1                        # Nodes to run on
#export PPN=4                          # Processors per node.
#export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
#export MEMUPPERBOUND=45000             # In KB


#export NPERSCRIPT=150    # Number of rows to run per script.
#export NNODES=2                        # Nodes to run on
#export PPN=8                          # Processors per node.
#export TIMEUPPERBOUND=1439              # In minutes. This is an upper bound for some of the elastic net fitting.
#export MEMUPPERBOUND=45000             # In KB


export SCRIPTNUM=$(expr $NROWS / $NPERSCRIPT + 1)
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