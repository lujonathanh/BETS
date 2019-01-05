#!/bin/bash


# Assume both of below files to be n genes x T timepoints
# RAND should be DATAFILE where each gene is randomized by time.

#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean/small_er-rep3-2cutoff-0mean-1var-100.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/small_data/0mean/small_er-rep3-2cutoff-0mean-1var-100-rand.txt
#export NROWS=100
#export GENES=test-new-fdr
#export LOADREPS=0



#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/small_data/0mean1var/small-100.CEBPB.rep_file_list.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/small_data/0mean1var/small-100.CEBPB.rep_file_list-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/small_data/filter_gene_list/gene_list.txt
#export NROWS=100
#export GENES=test.CEBPB..slow
#export NORMALIZATION=0mean1var
#export LOADREPS=1




#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/0mean1var/CEBPB.rep_file_list.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/0mean1var/CEBPB.rep_file_list-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/over-expression/CEBPB/filter_genelist/CEBPB.sig_genes_reg_fdr-0.05-all-ensg-cutoff-2-plusGR.txt
#export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is only used to determine how to parallelize
#export GENES=CEBPB           #protein-coding
#export NORMALIZATION=0mean1var
#export LOADREPS=1



#export DATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean1var/reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean1var/reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean1var/genes.txt
#export NROWS=100     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=            #protein-coding
#export NORMALIZATION=0mean1var
#export LOADREPS=1


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

#export DATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/DREAM/data/insilico_size100_1/0mean/genes.txt
#export NROWS=100     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=TEST-PAIRWISE-insilico_size100_1_bootstrap            #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1

#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/edgeR-reg-cutoff-2-plus-GR-0mean_joint-unperturbed-reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/edgeR-reg-cutoff-2-plus-GR-0mean_joint-unperturbed-reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/sig_genes_reg_fdr-0.05-all-ensg-cutoff-2-plusGR.txt
#export NROWS=2768     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=joint-unperturbed           #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1


#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/edgeR-reg-cutoff-2-plus-GR-0mean_joint-unperturbed-reps-urand.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/edgeR-reg-cutoff-2-plus-GR-0mean_joint-unperturbed-reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/unperturbed/data/joint-unperturbed/0mean/sig_genes_reg_fdr-0.05-all-ensg-cutoff-2-plusGR.txt
#export NROWS=2768     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=joint-unperturbed_urand           #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1


##
#export DATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean/reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean/reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/DREAM/data//0mean/genes.txt
#export NROWS=100     # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=TEST-LITE-BOOTSTRAP-METHOD-            #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1

#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/0mean1var/integration.rep_file_list.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/0mean1var/integration.rep_file_list-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/sig_genes_reg_fdr-0.05-all-ensg-cutoff-2-plusGR.txt
#export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=integration            #protein-coding
#export NORMALIZATION=0mean1var
#export LOADREPS=1


#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/0mean/edgeR-reg-cutoff-2-plus-GR-0mean-reps-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/raw_files/sig_genes_reg_fdr-0.05-all-ensg-cutoff-2-plusGR.txt
#export NROWS=2768      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=prot2TPM.bootstrap          #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1

#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/integration.rep_file_list.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/integration.rep_file_list-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/filter_gene_list/gene_list.txt
#export NROWS=100      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=testingsmallintegration            #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1


#export DATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/integration.rep_file_list.txt
#export RANDDATAFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/integration.rep_file_list-rand.txt
#export GENEFILE=/Users/jlu96/v-causal-snps/data/GeneExpressionData/edgeR-reg-reps-norm/data/final-integration/small_data/filter_gene_list/gene_list.txt
#export NROWS=100      # Rows (i.e. # genes) in your datafile. This is used to determine how to parallelize
#export GENES=testingsmallintegration            #protein-coding
#export NORMALIZATION=0mean
#export LOADREPS=1


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
#export CAUSAL=ridge
#export CAUSAL=lasso
#export LAG=2
export LAG=1
export NULL=g

export BOOTSTRAPNUM=1000
#export BOOTSTRAPNUM=100
#export BOOTSTRAPNUM=10
# you have to run prep_jobs_bootstrap.sh in order to get the bootstrap jobs



# Output names
export OUTPUTNAME="$GENES-$NORMALIZATION-$REPS-$DEG-$SAMPLE-$CAUSAL-$LAG-$NULL"
export FOLDER=della/$OUTPUTNAME                        # Folder made on your directory to scp to cluster.


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
